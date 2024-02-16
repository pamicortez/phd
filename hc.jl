function addNode(k, t, u, prio, dday, availableTime, sortedDPprio, visited)
    j = ""
    if u == o
        i = 1
        if FIRSTCHOICE == BYPRIORITY
            while i <= nDP && ((sortedDPprio[i] in visited) || (distance[o,sortedDPprio[i],t]+distance[sortedDPprio[i],o,t]) + vloadtime[vmodel[k]] + unloadtime >
                vspeed[vmodel[k]]*availableTime)
                i += 1
            end
            if i <= nDP && !(sortedDPprio[i] in visited) && (distance[o,sortedDPprio[i],t]+distance[sortedDPprio[i],o,t]) + vloadtime[vmodel[k]] + unloadtime <=
                vspeed[vmodel[k]]*availableTime
                j = sortedDPprio[i]
            end
        else # BYDEMAND
            biggest = 0
            for jaux in DP
                d = sum(priority[c] * dday[c,jaux] for c in C)
                if biggest < d && !(jaux in visited) && (distance[o,jaux,t]+distance[jaux,o,t]) + vloadtime[vmodel[k]] + unloadtime <= vspeed[vmodel[k]]*availableTime
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
    return j
end

# transfer cargo to smaller truck, and if possible, an already rented one
function getBestTruck(t, originalk, route, routetime, routeweight, availableTime, rent)
    bestk = bestt = -1
    bestvmodel = getBestVmodel(t, t, originalk, route, routeweight, availableTime)
    if bestvmodel == ""
        bestvmodel = vmodel[originalk]
    end
    # println(bestvmodel, " originalk ", originalk, " ", t1, " ", t2)

    # for t in t1:t2
        for k2 in nK:-1:1
            k = K[k2]
            if (vmodel[k] == bestvmodel && routetime-ZERO <= availableTime[k,t]) &&
                (bestk == -1 || (rent[k,t] > 1 && availableTime[k,t] < availableTime[bestk,t]))
                bestk = k
                # bestt = t
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

# remove y from route | routeloads[y] = minimum(routeloads) and y not in [o,i,j]
function remove_element(i, j, k, t, p, route, qdict, available_time)
    routeloads = get_route_load(k,t,p, route, qdict)
    delete!(routeloads, o)
    delete!(routeloads, i)
    delete!(routeloads, j)

    while available_time < getRouteDistance(route,t)/vspeed[vmodel[k]] + (length(route)-2)*unloadtime + vloadtime[vmodel[k]]
        minrouteload = minimum(routeloads)
        y = findall(x->x==minrouteload, routeloads)[1]
        index = findall(x->x==y, route)[1]
        deleteat!(route, index)
        delete!(routeloads, y)
    end
end

function get_maximum_load(k,t,p, newnodes, originalroute, route, qdict, loads)
    fo = 0
    available = get_available_load(k,t,p, originalroute, route, qdict)
    unsatisfied = Dict(c => sum(demand[c,y,tau] - loads[c,y,tau] for y in newnodes
         for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t) for c in C)

    if sum(values(unsatisfied)) > 1
        for c in C
            if available < 1
                break
            end
            aux = min(unsatisfied[c], available/weight[c])
            if aux > 1
                fo -= priority[c] * aux
                available -= weight[c] * aux
            end
        end
    end

    return fo
end

function add_arc(route, i,j,t)
    mindist = B
    minindex = 0
    for index in 1:length(route)-1
        d = distance[route[index],i,t] + distance[j,route[index+1],t] - distance[route[index],route[index+1],t]
        if mindist > d
            mindist = d
            minindex = index+1
        end
    end
    if minindex > 0
        insert!(route, minindex, i)
        insert!(route, minindex+1, j)
    end
end

function delete_route(k,t,p, vset, qdict, fuel, loads)
    i = o
    j = ""
    for j2 in DP
        if (i,j2,k,t,p) in vset
            j = j2
            break
        end
    end
    while j != o
        for c in C
            loads[c,j,t] -= qdict[c,i,j,k,t,p]
            qdict[c,i,j,k,t,p] = 0
        end
        delete!(vset, (i,j,k,t,p))
        fuel[t] -= vfuel[vmodel[k]] * distance[i,j,t]

        j2 = ""
        for j3 in N
            if (j,j3,k,t,p) in vset
                j2 = j3
                for c in C
                    loads[c,j,t] += get(qdict, (c,j,j2,k,t,p), 0)
                end
                break
            end
        end

        i = j
        j = j2
    end
    delete!(vset, (i,j,k,t,p))
    fuel[t] -= vfuel[vmodel[k]] * distance[i,j,t]
end

function destroy_random_routes(vset, qdict, rent, fuel, loads)
    lastp = Dict()
    if rand() > 0.5
        for t in T
            for k in K
                if ("theta",o,k,t,1) in vset # vehicle k is in use day t
                    pmax = 0
                    for p in 3:nP
                        if (o,"theta",k,t,p) in vset
                            delete!(vset, (o,"theta",k,t,p))
                            if p == 3
                                delete!(vset, ("theta",o,k,t,1))
                                rent[k,t] = 0
                            else
                                push!(vset, (o,"theta",k,t,p-1))
                            end
                            pmax = p-1
                            break
                        end
                    end
                    lastp[k,t] = p = pmax
                    delete_route(k,t,p, vset, qdict, fuel, loads)
                end
            end
        end
    else
        t = 0
        k = ""
        while !(("theta",o,k,t,1) in vset)
            t = rand(1:nT)
            k = "k"*string(rand(1:nK))
        end
        rent[k,t] = 0
        lastp[k,t] = 0
        for p in 2:nP
            if !((o,"theta",k,t,p) in vset)
                delete_route(k,t,p, vset, qdict, fuel, loads)
            else
                delete!(vset, ("theta",o,k,t,1))
                delete!(vset, (o,"theta",k,t,p))
                break
            end
        end
    end
    return lastp
end

function destroy_create(vset, qdict, rent, fuel, loads, incumbent)
    start = time()
    elapsedtime = 0
    psfo = psrent = psfuel = B
    v = q = unsatisfied_demand = nothing

    while elapsedtime < MY_TIME_LIMIT*10
        lastp = destroy_random_routes(vset, qdict, rent, fuel, loads)

        available_time = Dict()
        for t in T
            for k in K
                available_time[k,t] = get_available_time(k, t, vset)
            end
        end

        # create up to T.K new routes
        for k in K
            for t in T
                if (k,t) in keys(lastp)
                    route = ["a","b","c"]
                    while available_time[k,t] > vloadtime[vmodel[k]] && length(route) > 1
                        residualsupply, residualdemand = getResidual(nT+1, loads)
                        route = saveIncumbent(k, t, residualsupply, residualdemand, vset, qdict, rent, fuel, loads, lastp, available_time)
                    end
                end
            end
        end

        residualsupply, residualdemand = getResidual(nT+1, loads)
        auxfo   = sum(priority[c] * residualdemand[c,j,t] for c in C for j in DP for t in T)
        auxrent = sum(values(rent))
        auxfuel = sum(values(fuel))
        if incumbent > auxfo + auxrent + auxfuel
            psfo, psrent, psfuel = auxfo, auxrent, auxfuel
            v = copy(vset)
            q = copy(qdict)
            unsatisfied_demand = copy(residualdemand)
        end

        elapsedtime = time() - start
    end

    if DEBUG && v !== nothing
        printRouting(C, No, K, T, P, v, q, vload, vmodel, unsatisfied_demand, "DC", psfo + psrent + psfuel) # +eq
        flush(stdout)
    end

    return psfo, psrent, psfuel
end
