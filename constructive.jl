function createSolution(method, t1, t2, demandhc, supplyhc)
    vset, qdict = Set(), Dict()
    loads = Dict()
    rent, fuel = Dict{Any, Float64}(), Dict{Integer, Float64}()
    availableTime, lastp = Dict{Any, Float64}(), Dict{Any, Integer}()

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

    elapsedtime = 0
    if method != "notrip"
        elapsedtime = @elapsed fullLoadTrucks(t1:t2, supplyhc, demandhc, vset, qdict, rent, fuel, loads, availableTime, lastp, k1, k2, ordering)
    else
        setVars(loads, fuel, rent, availableTime, lastp)
    end
    elapsedtime += @elapsed routing(t1, t2, supplyhc, demandhc, vset, qdict, rent, fuel, loads, availableTime, lastp)
    eq = getEquity(demandhc)

    return demandhc, loads, rent, fuel, sum(values(eq)), Statistics.mean(values(eq)), maximum(values(eq)), vset, qdict, elapsedtime
end

function fullLoadTrucks(T, supply, demand, vset, qdict, rent, fuel, loads, availableTime, lastp, k1, k2, ordering)
    setVars(loads, fuel, rent, availableTime, lastp)
    aux = Dict()
    lastnode = 0

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
                    roundTrip = vloadtime[vmodel[k]] + (distance[o,j,t] + distance[j,o,t])/vspeed[vmodel[k]]
                    if availableTime[k,t] >= roundTrip

                        toDeliver = 0
                        for c in C
                            toDeliver += weight[c] * min(sum(supply[c,t2] for t2 in 1:t),
                                sum(demand[c,j,t2] for t2 in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t))
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

                availableTime[k,t] -= vloadtime[vmodel[k]] + (distance[o,j,t] + distance[j,o,t])/vspeed[vmodel[k]]
                push!(vset, (o,j,k,t,p))
                push!(vset, (j,o,k,t,p))
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

function routing(t1, t2, supply1, demand1, vset, qdict, rent, fuel, loads, availableTime, lastp)
    t = t1
    while t <= t2
        improved = true
        while improved
            improved = visitingByPrio(t, tend[C[1],t], supply1, demand1, vset, qdict, rent, fuel, loads, availableTime, lastp)
        end
        t = tend[C[1],t]+1
    end
end

function visitingByPrio(t1, t2, supply1, demand1, vset, qdict, rent, fuel, loads, availableTime, lastp)
    bestfo = bestk = bestt = B

    for t in t1:t2
        sday, dday = preprocessing(t, supply1, demand1)
        prio = getPrio(t, sday, dday, loads)
        if sum(values(prio)) > ZERO

            sortedDPprio = getDPsortedByPrio(prio, dday)
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
            # println(t, " big truck ", k)
            if k == "" # no vehicle available at time t
                continue
            end

            route, q, routeload, routetime, routeweight = createRoute(k, t, prio, sday, dday, availableTime[k,t], sortedDPprio)
            if routeweight > 1

                k = getBestTruck(t, k, route, routetime, routeweight, availableTime, rent)
                # println(route, " ", t, " best truck ", k)

                fo = 0
                if rent[k, t] < ZERO
                    fo += vcost[vmodel[k]]
                end
                u = o
                for j in route
                    fo += vfuel[vmodel[k]] * distance[u,j,t]
                    if j != o
                        itens = 0
                        for c in C
                            fo -= prio[c,j] * q[c,j,t]
                            itens += 1
                        end
                        # fo -= itens * sum(q[c,u,j] for c in C)/10 # good solutions deliver 2+ relief itens
                    end
                    u = j
                end
                if bestfo-ZERO > fo
                    bestfo = fo
                    bestk = k
                    bestt = t
                end
            end
        end
    end

    if bestt <= nT
        # println(bestt, " saveIncumbent ", bestk)
        saveIncumbent(bestk, bestt, supply1, demand1, vset, qdict, rent, fuel, loads, lastp, availableTime)
    end

    return bestt <= nT
end

function createRoute(k, t, prio, sday, dday, availableTime, sortedDPprio)
    visited, route = Array{String}(undef, 0), Array{String}(undef, 0)
    routetime = vloadtime[vmodel[k]]
    routeweight = 0
    q, routeload = Dict(), Dict()
    for c in C
        routeload[c] = 0
    end

    u = o
    routeweight, u = addNodes(k, t, prio, sday, dday, availableTime, u, visited, route, routetime, routeweight, q, routeload, sortedDPprio)

    if u != o
        push!(route, o)
        if length(route) > 3 && length(route) < 10
            elapsedtime = @elapsed route = tsp(route, t)
            # println(t, " tsp", route, " ", elapsedtime)
        end
        # if length(route) >= 10 println(t, " before ", route) end
        routetime = vloadtime[vmodel[k]] + getRouteDistance(route, t)/vspeed[vmodel[k]]
    end

    return route, q, routeload, routetime, routeweight
end

function addNodes(k, t, prio, sday, dday, availableTime, u, visited, route, routetime, routeweight, q, routeload, sortedDPprio)
    while length(visited) < nDP && routeweight < vload[vmodel[k]]
        j = ""
        if u == o
            i = 1
            if FIRSTCHOICE == BYPRIORITY
                while i <= nDP && ((sortedDPprio[i] in visited) || (distance[o,sortedDPprio[i],t]+distance[sortedDPprio[i],o,t]) > vspeed[vmodel[k]]*availableTime)
                    i += 1
                end
                if i <= nDP && !(sortedDPprio[i] in visited) && (distance[o,sortedDPprio[i],t]+distance[sortedDPprio[i],o,t]) <= vspeed[vmodel[k]]*availableTime
                    j = sortedDPprio[i]
                end
            else # BYDEMAND
                biggest = 0
                for jaux in DP
                    d = sum(priority[c] * dday[c,jaux] for c in C)
                    if biggest < d && !(jaux in visited) && (distance[o,jaux,t]+distance[jaux,o,t]) <= vspeed[vmodel[k]]*availableTime
                        biggest = d
                        j = jaux
                    end
                end
            end
        else
            sortedDP = nearestNeighbours(t, u, prio)
            # println(u, " nearestNeighbours ", sortedDP)
            i = 1
            while i <= nDP && (sortedDP[i] in visited || sum(dday[c,sortedDP[i]] for c in C) < ZERO)
                i += 1
            end
            if i <= nDP && !(sortedDP[i] in visited) && sum(dday[c,sortedDP[i]] for c in C) > ZERO
                j = sortedDP[i]
            end
            # println(u, " - ", j, " - ", visited, " ", j in visited)
        end
        if j == ""
            break
        end
        push!(visited, j)

        if routetime + (distance[u,j,t] + distance[j,o,t])/vspeed[vmodel[k]] < availableTime

            # verificar se load de j pode ser entregue em única entrega, mas estaria sendo "quebrado" em 2 entregas
            currentload = 0
            for c in C
                currentload += min(dday[c,j], sday[c]) * weight[c]
            end
            if currentload % vload[vmodel[k]] < vload[vmodel[k]] - routeweight + 0.1

                for c in C
                    currentload = min(dday[c,j], sday[c], (vload[vmodel[k]] - routeweight)/weight[c])
                    # println(j, " ", c, " ", dday[c,j], " ", sday[c], " ", currentload)
                    if currentload > 1
                        q[c,j,t] = currentload
                        routeload[c] += currentload
                        routeweight += currentload*weight[c]
                    else
                        q[c,j,t] = 0
                    end
                end

                if sum(q[c,j,t] for c in C) > 1
                    routetime += distance[u,j,t] / vspeed[vmodel[k]]
                    push!(route, j)
                    u = j
                end
            end
        end
    end

    return routeweight, u
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

    return getBestTruck(t1, t, k, route, routeweight, availableTime, rent)
end

# transfer cargo to smaller truck, and if possible, an already rented one
function getBestTruck(t, originalk, route, routetime, routeweight, availableTime, rent)
    bestk = bestt = -1
    bestvmodel = getBestVmodel(t, t, originalk, route, routeweight, availableTime)
    if bestvmodel == ""
        bestvmodel = vmodel[originalk]
    end
    bestrented = false
    # println(bestvmodel, " originalk ", originalk, " ", t1, " ", t2)

    # for t in t1:t2
        for k2 in nK:-1:1
            k = K[k2]
            if vmodel[k] == bestvmodel && routetime-ZERO <= availableTime[k,t]
                if bestt == -1 || rent[k, t] > 1
                    if !bestrented || (rent[k, t] > 1 && availableTime[k,t] < availableTime[bestk,t])
                        bestk = k
                        # bestt = t
                    end
                end
            end
        end
    # end

    if bestk == -1
        # for t in t1:t2
        #     if validRouteTime(route, originalk, t, availableTime)
        #         return originalk, t
        #     end
        # end
        return originalk
    end
    return bestk#, bestt
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

    return vloadtime[vmodel[k]] + dist/vspeed[vmodel[k]] <= availableTime[k,t]
end

function getBigTruckWithMaximumAvailableTime(t, j, availableTime)
    maxAvailableTime = 0
    model = k = ""

    for kaux in K # K is ordered by vload
        if (model == "" || vmodel[kaux] == model) && vehicle[kaux,t] && availableTime[kaux,t] >=
            vloadtime[vmodel[kaux]] + (distance[o,j,t] + distance[j,o,t])/vspeed[vmodel[kaux]]
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
    route, q, routeload, routetime, routeweight = createRoute(k, t, prio, sday, dday, availableTime[k,t], getDPsortedByPrio(prio,dday))

    rent[k,t] = vcost[vmodel[k]]
    if lastp[k,t] == 0
        lastp[k,t] = 2
        push!(vset, ("theta",o,k,t,1))
    else
        delete!(vset, (o,"theta",k,t,lastp[k,t]))
    end
    push!(vset, (o,"theta",k,t,lastp[k,t]+1))

    u = o
    for j in route
        push!(vset, (u,j,k,t,lastp[k,t]))
        fuel[t] += vfuel[vmodel[k]] * distance[u,j,t]

        for c in C
            if j != o
                qdict[c,u,j,k,t,lastp[k,t]] = q[c,j,t] # deltaQ
                loads[c,j,t] += q[c,j,t]
            end
        end

        if j != o
            updateResidual(C, [j], t, supply1, demand1, q)
            u = j
        end
    end

    # println(availableTime[k,t], " ", route)
    availableTime[k,t] -= routetime
    # println(availableTime[k,t], " ", routetime)
    lastp[k,t] += 1
end
