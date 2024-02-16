function createSolution(method, t1, t2, demandhc, supplyhc)
    vset, qdict, loads = Set(), Dict(), Dict()
    rent, fuel = Dict{Any, Float64}(), Dict{Integer, Float64}()
    availableTime, lastp = Dict{Any, Float64}(), Dict{Any, Integer}()
    setVars(loads, fuel, rent, availableTime, lastp)

    routes = Dict()
    for k in K
        for t in t1:t2
            routes[k,t] = []
        end
    end

    if method != "notrip"
        fullLoadTrucks(t1:t2, supplyhc, demandhc, vset, qdict, rent, fuel, loads, availableTime, lastp, method, routes)
    end

    routing(t1, t2, supplyhc, demandhc, vset, qdict, rent, fuel, loads, availableTime, lastp, routes)
    eq = getEquity(demandhc)

    return routes, demandhc, rent, fuel, sum(values(eq)), vset, qdict, loads
end

function getOrdering(method)
    k1, k2, ordering = 0, 0, ""
    if method == "1:nK"
        k1, k2, ordering = 1, nK, 1
    elseif method == "nK:1"
        k1, k2, ordering = nK, 1, -1
    elseif method == "onlybig"
        k1, k2, ordering = 1, ceil(Int, nK/2), 1
    elseif method == "onlysmall"
        k1, k2, ordering = nK, ceil(Int, nK/2)+1, -1
    elseif method == "both"
        k1, k2, ordering = ceil(Int, nK/4)+1, floor(Int, 3*nK/4), 1
    end
    return k1, k2, ordering
end

function fullLoadTrucks(T, supply, demand, vset, qdict, rent, fuel, loads, availableTime, lastp, method, routes)
    aux = Dict()
    lastnode = 0

    k1, k2, ordering = getOrdering(method)
    for kaux in k1:ordering:k2
        k = K[kaux]
        jj = firstj = lastnode == nDP-1 ? nDP : mod(lastnode+1, nDP)

        improvement = false
        while true
            j = DP[jj]

            bestfo = bestt = B
            empty!(aux)
            for t in T
                if vehicle[k,t]
                    roundTrip = vloadtime[vmodel[k]] + unloadtime + (distance[o,j,t] + distance[j,o,t])/vspeed[vmodel[k]]
                    if availableTime[k,t] >= roundTrip

                        toDeliver = 0
                        for c in C
                            toDeliver += weight[c] * min(sum(supply[c,t2] for t2 in 1:t), sum(demand[c,j,t2] for t2 in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t))
                        end
                        if toDeliver >= vload[vmodel[k]]
                            delivered = 0
                            fo = vfuel[vmodel[k]] * (distance[o,j,t] + distance[j,o,t])
                            for c in C # load vehicle
                                if vload[vmodel[k]]-delivered > weight[c]
                                    aux[c,j,t] = min((vload[vmodel[k]]-delivered)/weight[c], sum(supply[c,t2] for t2 in 1:t),
                                        sum(demand[c,j,t2] for t2 in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t))
                                    delivered += weight[c] * aux[c,j,t]
                                    fo -= priority[c] * aux[c,j,t]
                                else
                                    aux[c,j,t] = 0
                                end
                            end

                            if bestfo-ZERO > fo
                                bestfo = fo
                                bestt = t
                            end
                        end
                    end
                end
            end

            if bestt <= nT
                improvement = true
                lastnode = jj
                t = bestt

                if lastp[k,t] == 0
                    push!(vset, ("theta",o,k,t,1))
                    push!(vset, (o,"theta",k,t,3))
                    rent[k,t] = vcost[vmodel[k]]
                    lastp[k,t] = 3
                else
                    delete!(vset, (o,"theta",k,t,lastp[k,t]))
                    lastp[k,t] += 1
                    push!(vset, (o,"theta",k,t,lastp[k,t]))
                end
                p = lastp[k,t]-1

                for c in C
                    qdict[c,o,j,k,t,p] = aux[c,j,t]
                    loads[c,j,t] += aux[c,j,t]
                    updateResidual([c], [j], t, supply, demand, aux)
                end
                # c=C[1]
                # println(bestt, " ", sum(supply[c,t2] for t2 in 1:t), " ", sum(demand[c,j,t2] for j in DP for t2 in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t))

                availableTime[k,t] -= vloadtime[vmodel[k]] + unloadtime + (distance[o,j,t] + distance[j,o,t])/vspeed[vmodel[k]]
                push!(vset, (o,j,k,t,p))
                push!(vset, (j,o,k,t,p))
                push!(routes[k, t], [o,j,o])
                fuel[t] += vfuel[vmodel[k]] * (distance[o,j,t] + distance[j,o,t])
            end

            jj = jj == nDP-1 ? nDP : mod(jj+1, nDP)
            if jj == firstj
                if !improvement
                    break
                end
                improvement = false
            end
        end
    end
end

function routing(t1, t2, supply1, demand1, vset, qdict, rent, fuel, loads, availableTime, lastp, routes)
    tm = t1
    while tm <= t2
        while true
            bestfo = bestk = bestt = B-10
            for t in tm:tend[C[1],tm]
                fo, k, route = createTentativeRoute(t, supply1, demand1, vset, qdict, rent, fuel, loads, availableTime, lastp)
                if bestfo > fo
                    bestfo, bestk, bestt = fo, k, t
                    # println(fo, " ", k, " ", t, " ", route)
                end
            end
            if bestt > nT
                break
            end
            route = saveIncumbent(bestk, bestt, supply1, demand1, vset, qdict, rent, fuel, loads, lastp, availableTime)
            push!(routes[bestk, bestt], route)
            # println(bestk, " ", bestt, " route ch ", route)
            # println(routes[bestk, bestt])
        end
        tm = tend[C[1],tm]+1
    end
end

function createTentativeRoute(t, supply1, demand1, vset, qdict, rent, fuel, loads, availableTime, lastp)
    fo, k, r = B, "", []
    sday, dday = preprocessing(t, supply1, demand1)
    if sum(values(dday)) > ZERO
        prio = getPrio(t, sday, dday, loads)
        sortedDPprio = getDPsortedByPrio(prio, dday)
        k = getBigTruck(t, sortedDPprio, availableTime)
        # println(t, " big truck ", k)

        if k != "" # else: no vehicle available at time t
            route, q, routeload, routeweight, routetime, routedistance = createRoute(k, t, prio, sday, dday, availableTime[k,t], sortedDPprio)
            if routeweight > 1
                k = METHOD == 6 ? getBestTruckG(t, k, route, routetime, routeweight, availableTime, rent) : getBestTruck(t, k, route, routetime, routeweight, availableTime, rent)
                # println(route, " ", t, " best truck ", k)
                fo = (rent[k,t] < 1 ? vcost[vmodel[k]] : 0) + vfuel[vmodel[k]] * routedistance -
                    sum(prio[c,j] * q[c,j,t] for c in C for j in filter(e->e!=o, route))
                r = route
            end
        end
    end
    return fo, k, r
end

function createRoute(k, t, prio, sday, dday, availableTime, sortedDPprio)
    visited, route = Array{String}(undef, 0), Array{String}(undef, 0)
    routetime = vloadtime[vmodel[k]]
    q, routeload = Dict(), Dict()
    for c in C
        routeload[c] = 0
    end
    routedistance = routeweight = 0
    push!(route, o)
    u = o

    while length(visited) < nDP && routeweight < vload[vmodel[k]]
        j = METHOD == 6 ? addNodeG(k, t, u, prio, dday, availableTime, sortedDPprio, visited) : addNode(k, t, u, prio, dday, availableTime, sortedDPprio, visited)
        if j == ""
            break
        end
        push!(visited, j)

        if routetime + unloadtime + (distance[u,j,t] + distance[j,o,t])/vspeed[vmodel[k]] < availableTime
            # verificar se load de j pode ser entregue em única entrega, mas estaria sendo "quebrado" em 2 entregas
            currentload = sum(min(dday[c,j], sday[c]) * weight[c] for c in C)
            if currentload % vload[vmodel[k]] < vload[vmodel[k]] - routeweight + 0.1
                for c in C
                    currentload = min(dday[c,j],sday[c],(vload[vmodel[k]] - routeweight)/weight[c])
                    # println(j, " ", c, " ", dday[c,j], " ", sday[c], " ", currentload)
                    if currentload > 1
                        q[c,j,t] = currentload
                        routeload[c] += currentload
                        routeweight += currentload*weight[c]
                        sday[c] -= currentload
                    else
                        q[c,j,t] = 0
                    end
                end

                if sum(q[c,j,t] for c in C) > 1
                    routetime += distance[u,j,t] / vspeed[vmodel[k]] + unloadtime
                    push!(route, j)
                    u = j
                end
            end
        end
    end

    if u != o
        push!(route, o)
        if length(route) > 4
            route = length(route) < 10 ? tsp(route, t) : getOrderByNearest(route, t)
        end
        routedistance = getRouteDistance(route,t)
        routetime = vloadtime[vmodel[k]] + (length(route)-2)*unloadtime + routedistance/vspeed[vmodel[k]]
    end

    return route, q, routeload, routeweight, routetime, routedistance
end

function scheduleDelivery(t, k, supply1, demand1, route, routeload, routeweight, availableTime, rent)
    # bring forward the delivery -- no! it's better to test for all t, because road conditions change everyday
    t1 = t
    auxt = 1
    for c in C
        if routeload[c] > 0
            auxt = max(auxt, div(t-1,threshold[c]+1)*(threshold[c]+1)+1)
        end
    end
    for taux in auxt:t-1
        valid = true
        sday, dday = preprocessing(taux, supply1, demand1)

        for c in C
            if sday[c] < routeload[c]-ZERO || sum(dday[c,j] for j in DP) < routeload[c]-ZERO
                valid = false
                break
            end
        end

        if valid
            t1 = taux
            break
        end
    end

    return METHOD == 6 ? getBestTruckG(t, k, route, routetime, routeweight, availableTime, rent) : getBestTruck(t1, t, k, route, routeweight, availableTime, rent)
end

function getBestVmodel(t1, t2, originalk, route, routeweight, availableTime)
    for k2 in nK:-1:findall(x->x==originalk, K)[1]+1
        k = K[k2]
        if vload[vmodel[k]] >= routeweight-0.1
            for t in t1:t2
                if validRouteTime(route, k, t, availableTime)
                    return vmodel[k]
                end
            end
        end
    end
    return ""
end

function validRouteTime(route, k, t, availableTime)
    if !vehicle[k,t]
        return false
    end

    dist = 0
    u = o
    for j in route
        dist += distance[u,j,t]
        u = j
    end
    dist += distance[u,o,t]

    return vloadtime[vmodel[k]] + length(route)*unloadtime + dist/vspeed[vmodel[k]] <= availableTime[k,t]
end

function getBigTruck(t, sortedDPprio, availableTime)
    i = 1
    k = ""
    while k == "" && i <= nDP
        while i <= nDP && distance[o,sortedDPprio[i],t] > 400
            # println(i, " ", sortedDPprio[i])
            i += 1
        end
        if i <= nDP
            k = getBigTruckWithMaximumAvailableTime(t, sortedDPprio[i], availableTime)
            # println(t, " ", k, " ", sortedDPprio[i])
        end
        i += 1
    end
    return k
end

function getBigTruckWithMaximumAvailableTime(t, j, availableTime)
    maxAvailableTime = 0
    model = k = ""

    for kaux in K # K is ordered by vload
        if (model == "" || vmodel[kaux] == model) && vehicle[kaux,t] && availableTime[kaux,t] >=
            vloadtime[vmodel[kaux]] + unloadtime + (distance[o,j,t] + distance[j,o,t])/vspeed[vmodel[kaux]]
            if maxAvailableTime < availableTime[kaux,t]
                maxAvailableTime = availableTime[kaux,t]
                k = kaux
                model = vmodel[kaux]
            end
        end
    end

    return k
end

function saveIncumbent(k, t, supply1, demand1, vset, qdict, rent, fuel, loads, lastp, availableTime)
    sday, dday = preprocessing(t, supply1, demand1)
    prio = getPrio(t, sday, dday, loads)
    route, q, routeload, routeweight, routetime, routedistance = createRoute(k, t, prio, sday, dday, availableTime[k,t], getDPsortedByPrio(prio,dday))

    if length(route) > 2
        rent[k,t] = vcost[vmodel[k]]
        if lastp[k,t] == 0
            lastp[k,t] = 2
            push!(vset, ("theta",o,k,t,1))
        elseif (o,"theta",k,t,lastp[k,t]) in vset
            delete!(vset, (o,"theta",k,t,lastp[k,t]))
        end
        closed = false
        for p in lastp[k,t]+1:nP
            if (o,"theta",k,t,p) in vset
                closed = true
                break
            end
        end
        p = lastp[k,t]
        if !closed
            push!(vset, (o,"theta",k,t,p+1))
        end

        u = o
        for j in route[2:end]
            push!(vset, (u,j,k,t,p))
            fuel[t] += vfuel[vmodel[k]] * distance[u,j,t]

            if j != o
                for c in C
                    # println("q ", u, " ", j, " ", p)
                    qdict[c,u,j,k,t,p] = q[c,j,t] # deltaQ
                    loads[c,j,t] += q[c,j,t]
                end

                updateResidual(C, [j], t, supply1, demand1, q)
                u = j
            end
        end

        if length(route) >= 4
            # println(route)
            u = route[end-1]
            j = route[end-2]
            i = ""
            for i in reverse(route[1:end-3])
                # println(i, " ", j, " ", u, " ", p)
                for c in C
                    qdict[c,i,j,k,t,p] += qdict[c,j,u,k,t,p]
                end
                u = j
                j = i
            end
        end

        # println(availableTime[k,t], " ", route)
        availableTime[k,t] -= routetime
        # println(availableTime[k,t], " ", routetime)
        lastp[k,t] += 1
    end

    return route
end
