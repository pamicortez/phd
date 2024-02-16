# latitude and longitude in degrees
function angleFromCoordinate(lat1, long1, lat2,long2)
    dlong = long2 - long1
    return 90-rad2deg(atan(sind(dlong)*cosd(lat2), cosd(lat1)*sind(lat2) - sind(lat1)*cosd(lat2)*cosd(dlong)))
end

# return cities ordered by angle from depot
function getCitiesByAngle()
    olat, olong = coordinate[o]

    angle = Array{Float64,1}(undef, nDP)
    dict = Dict{Any,Float64}()
    i = 0
    for j in DP
        i += 1
        lat, long = coordinate[j]
        angle[i] = angleFromCoordinate(olat, olong, lat,long)
        dict[j] = angle[i]
        # println("\"", j, "\"@", angle[i])
    end

    return DP[sort([1:nDP;], by=i->angle[i])] #,angle
end

function attractionVertices(t)
    max = 0
    besti = bestj = bestk = ""
    for i in DP
        for j in filter(e -> e != i, DP)
            for k in setdiff(DP, [i,j])
                sum = distance[i,j,t] + distance[j,k,t] + distance[k,i,t] # distance[o,i] + distance[o,j] + distance[o,k] +
                # println(i, " ", j, " ", k, " ", sum, " ", max)
                if max < sum
                    max = sum
                    besti, bestj, bestk = i, j, k
                end
            end
        end
    end

    indexi = findfirst(isequal(besti), DP)
    indexj = findfirst(isequal(bestj), DP)
    indexk = findfirst(isequal(bestk), DP)
    if indexi < indexj && indexi < indexk
        return indexj < indexk ? (besti, bestj, bestk) : (besti, bestk, bestj)
    elseif indexj < indexi && indexj < indexk
        return indexi < indexk ? (bestj, besti, bestk) : (bestj, bestk, besti)
    else
        return indexi < indexj ? (bestk, besti, bestj) : (bestk, bestj, besti)
    end
end

function attractionSets(i, j, k, t)
    Aij, Ajk, Aki = [], [], []
    push!(Aij, i)
    push!(Aij, j)
    push!(Ajk, j)
    push!(Ajk, k)
    push!(Aki, k)
    push!(Aki, i)

    for l in setdiff(N, [i,j,k])
        if distance[l,i, t] < distance[l,k, t] && distance[l,j, t] < distance[l,k, t]
            push!(Aij, l)
        elseif distance[l,j, t] < distance[l,i, t] && distance[l,k, t] < distance[l,i, t] # distance[l,i] < distance[l,o] && distance[l,j] < distance[l,o]
            push!(Ajk, l)
        else
             push!(Aki, l)
        end
    end

    # println(Aij, " ", Ajk, " ", Aki)
    return Aij, Ajk, Aki
end

# Cheapest Insertion: Find the city which insertion in the tour causes the smallest increase in length,
# i.e the city k which minimizes d(i, k)  + d(k, j) - d(i, j) with (i, j) an edge in the partial tour.
# https://github.com/afourmy/pyTSP
function cheapestInsertion(A, t)
    route = [A[1], A[2]]
    popfirst!(A)
    popfirst!(A)
    n = length(A)

    for i in 1:n
        min = B
        indexk = indexj = -1
        for j in 2:length(route)
            for k in 1:length(A)
                if min > distance[route[j-1],A[k], t] + distance[A[k],route[j], t] - distance[route[j-1],route[j], t]
                    min = distance[route[j-1],A[k], t] + distance[A[k],route[j], t] - distance[route[j-1],route[j], t]
                    indexk = k
                    indexj = j
                end
            end
        end
        insert!(route, indexj, A[indexk])
        deleteat!(A, indexk)
    end
    return route
end

function getCities(t)
    # sectors, orderedDPs = getSectors()
    #
    # ini = 1
    # for sector in sectors
        i, j, k = attractionVertices(t)
        Aij, Ajk, Aki = attractionSets(i, j, k, t) # DP[ini:sector-1]
    #     ini = sector
    # end
    route = [cheapestInsertion(Aij, t); cheapestInsertion(Ajk, t); cheapestInsertion(Aki, t)]
    indexi = findfirst(isequal(o), route)
    route = [route[indexi:end];route[1:indexi-1]]
    # println(route)

    indexi = 0
    for i in route
        indexi += 1
        for j in indexi+1:length(route)
            if i == route[j]
                deleteat!(route, j)
                break
            end
        end
    end
    return route[2:end] #move([route;o], i, j, k)
end

function furthest(toskip, orderedDPs, t) #sini:send
    maxdistance = 0.0
    maxi, maxj = "", ""
    candidates = setdiff(orderedDPs, toskip) #orderedDPs[sini:send]
    ncandidates = length(candidates)

    if ncandidates == 1
        i = maxi = candidates[1]
        for j in orderedDPs #[sini:send]
            if maxi != j && maxdistance < distance[i,j, t]
                maxdistance = distance[i,j, t]
                maxj = j
            end
        end
    elseif ncandidates > 1
        for i in 1:ncandidates
            for j in i+1:ncandidates
                if maxdistance < distance[candidates[i],candidates[j], t]
                    maxdistance = distance[candidates[i],candidates[j], t]
                    maxi, maxj = candidates[i], candidates[j]
                end
            end
        end
    end

    return maxi, maxj
end

function getIndex(j, index1, route1, orderedDPs)
    indexj = findfirst(isequal(j), orderedDPs)
    index1c = copy(index1)
    route1c = copy(route1)
    i = 1
    while i <= length(index1) && index1[i] < indexj
        i += 1
    end
    insert!(index1c, i, indexj)
    insert!(route1c, i+1, j)
    return index1c, route1c
end

function splitNearFarNodes(DP, t)
    avg = Statistics.mean(distance[o,j,t] for j in DP)
    near, far = [], []
    for j in DP
        push!(distance[o,j,t] <= avg ? near : far, j)
    end
    return near, far
end

function addRoute(set, routes, t)
    minspeed = minimum(vspeed)[2]
    maxloadtime = maximum(vloadtime)[2]
    orderedDPs = intersect(getCitiesByAngle(), set)

    i = 1
    while i < length(set)
        j = orderedDPs[i] #indexj
        route = [o,j]
        i += 1

        while i <= length(set) && getRouteDistance([route;orderedDPs[i];o], t)/minspeed + maxloadtime + unloadtime * length(route) < working_hours
            push!(route, orderedDPs[i])
            i += 1
        end
        if length(route) > 2
            push!(routes, length(route) < 10 ? tsp([route;o],t) : getOrderByNearest([route;o],t))
        end
    end
end

function onePetal(routes, t)
    near, far = splitNearFarNodes(DP, t)

    addRoute(near, routes, t)
    addRoute(far, routes, t)
    addRoute(DP, routes, t)
end

function twoPetals(routes, t)
    served = []
    orderedDPs = getCities(t) # one big route
    minspeed = minimum(vspeed)[2]
    maxloadtime = maximum(vloadtime)[2]

    # sectors = getSectors(orderedDPs, dictangle)
    # sini = 1
    # for send in sectors

        while length(served) < nDP
            i, j = furthest(served, orderedDPs, t) #, sini, send-1)
            if i == "" || j == ""
                break
            end
            # println(i, " 2-petal ", j)
            push!(served, i)
            push!(served, j)

            route1, route2 = [o,i,o], [o,j,o]
            index1 = [findfirst(isequal(i), orderedDPs)]
            index2 = [findfirst(isequal(j), orderedDPs)]

            for j in setdiff(orderedDPs, served)
                index1c, route1c = getIndex(j, index1, route1, orderedDPs)
                index2c, route2c = getIndex(j, index2, route2, orderedDPs)
                route1length = getRouteDistance(route1c, t)
                route2length = getRouteDistance(route2c, t)

                if min(route1length, route2length)/minspeed < working_hours - maxloadtime + (length(route1)-2)*unloadtime
                    if route1length < route2length
                        route1, index1 = route1c, index1c
                    else
                        route2, index2 = route2c, index2c
                    end
                    push!(served, j)
                end
            end

            # addRoute(routes, route1[1:end-1]) #petalLS(route1, dictangle))
            # addRoute(routes, route2[1:end-1]) #petalLS(route2, dictangle))
            push!(routes, length(route1) < 10 ? tsp(route1, t) : getOrderByNearest(route1,t))
            push!(routes, length(route2) < 10 ? tsp(route2, t) : getOrderByNearest(route2,t))
        end
    #     sini = send
    # end
end

function singleTrip(routes)
    for j in DP
        push!(routes, [o,j,o])
    end
end

function removeTrips(routes)
    visited = []
    for route in routes
        if length(route) == 3
            push!(visited, route[2])
        end
    end
    for j in visited
        deleteat!(routes, findall(x->x==[o,j,o], routes))
    end
end


function runPetal(routes, whichone, routesch=nothing, nt = 1, routesvalue=Dict(), qr=Dict())
    if length(routes) == 0
        aux = []
        singleTrip(aux)

        for t in T
            routes[t] = []
            onePetal(routes[t], t)
            twoPetals(routes[t], t)
            removeTrips(routes[t])
            for route in aux
                push!(routes[t], route)
            end

            nroutes = length(routes[t])
            for r in 1:nroutes
                routetime = getRouteDistance(routes[t][r],t)/maximum(values(vspeed)) + (length(routes[t][r])-2)*unloadtime + minimum(values(vloadtime))
                for p in 2:trunc(Int, working_hours/routetime)
                    push!(routes[t], routes[t][r])
                end
            end
        end
    end

    routeindex = Dict()
    if routesch !== nothing
        for t in T
            for k in K
                routeindex[k,t] = length(routes[t])+1
                # println("\n--------------- ", k, " ", t, " ", routeindex[k,t])
                for route in routesch[k,t]
                    push!(routes[t], route)
                    # println(t, " ", length(routes[t]), " ", route)
                end
            end
        end
    end
    incumbent, unsatisfiedpenalty, transportationcosts, eqpenalty, vset, qdict, loads = getRoutes(routes, routesch, routeindex, nt, routesvalue, qr)

    if DEBUG  && incumbent > 1
        printRouting(C, No, K, T, P, vset, qdict, vload, vmodel, nothing, whichone, incumbent)
    end
    GC.gc()
    sleep(SLEEP)
    GC.gc()

    return incumbent, unsatisfiedpenalty, transportationcosts, eqpenalty, vset, qdict, loads
end

# function getSectors()
#     orderedDPs, dictangle = getCitiesByAngle()
#
#     deltatheta = []
#     for i in 2:length(orderedDPs)
#         push!(deltatheta, abs(dictangle[orderedDPs[i]] - dictangle[orderedDPs[i-1]]))
#     end
#     deltatheta = mean(deltatheta) + std(deltatheta)
#
#     sectors = []
#     for i in 2:length(orderedDPs)
#         if abs(dictangle[orderedDPs[i]] - dictangle[orderedDPs[i-1]]) > deltatheta
#             push!(sectors, i)
#         end
#     end
#
#     push!(sectors, nDP+1)
#     # println("sectors ", sectors)
#     return sectors, orderedDPs
# end
#
# function move(route, i, j, k)
#     indexi = findfirst(isequal(i), route)
#     indexj = findfirst(isequal(j), route)
#     indexk = findfirst(isequal(k), route)
#
#     visited = Set()
#     while true
#         min = B
#         index = bestindex = 0
#         for index1 in indexi+1:indexj-1
#             for index2 in [2:indexi;indexj+1:indexk]
#                 newvalue = distance[route[index1-1], route[index1+1]] + distance[route[index2-1], route[index1]] + distance[route[index1], route[index2]]
#                 value = distance[route[index1-1], route[index1]] + distance[route[index1], route[index1+1]] + distance[route[index2-1], route[index2]]
#                 # println(route[index1-1], " ", route[index1], " ", route[index1+1], " ", route[index2-1], " ", route[index2], " ", value)
#                 # println(route[index1-1], " ", route[index1+1], " ", route[index2-1], " ", route[index1], " ", route[index2], " ", newvalue, "\n")
#                 if !(aux in visited) && min > newvalue && newvalue < value
#                     bestindex = index2
#                     index = index1
#                 end
#             end
#         end
#         if index > 0
#             push!(visited, index)
#             l = route[index]
#             deleteat!(route, index)
#             if bestindex > index
#                 bestindex -= 1
#             end
#             # print(l, " ", route[bestindex], " ")
#             insert!(route, bestindex, l)
#             println(getRouteDistance(route))
#             if bestindex < indexi
#                 indexi += 1
#             else
#                 indexj -= 1
#             end
#         else
#             break
#         end
#     end
#
#     return route
# end
#
# function getMinRoutes()
#     maxtype, maxcapacity = maximum(vload)
#     mintype, mincapacity = minimum(vload)
#     return min(nDP, nvehicles[maxtype,2] + max(0, ceil(Int, (sum(demand[C[1],j,1] for j in DP) - nvehicles[maxtype,2]*maxcapacity) / mincapacity)))
# end
#
# function geMinMaxDistance(route, index1, index2)
#     return (maximum(distance[o,j] for j in route[index1:index2-1])+minimum(distance[o,j] for j in route[index1:index2-1]))/2
# end
#
# function moveNode(route, j, indexj, i)
#     if distance[o,j] < distance[o,route[i]]
#         deleteat!(route, indexj)
#         insert!(route, i, j)
#     end
# end
#
# # up == route's first nodes. index is the position of the first "far" node
# function moveUp(route, index)
#     if index > 3
#         # calculating this measure, because I dont't want to change angle sorting for too little change
#         avgdist = geMinMaxDistance(route, 2, index)
#
#         indexj = index
#         for j in reverse(route[2:index-1])
#             indexj -= 1
#
#             i = indexj-1
#             while i > 1 && distance[o,j] < avgdist && distance[o,route[i]] > avgdist
#                 println(avgdist, " ", j, " ", distance[o,j], " ", route[i], " ", distance[o,route[i]])
#                 i -= 1
#             end
#             i += 1
#             # println(j, " ", i, " ", route[i])
#             moveNode(route, j, indexj, i)
#         end
#     end
# end
#
# function moveDown(route, index)
#     if index < length(route) - 2
#         avgdist = geMinMaxDistance(route, index+1, length(route))
#
#         indexj = index
#         for j in route[index+1:length(route)-1]
#             indexj += 1
#
#             i = indexj+1
#             while i < length(route) && distance[o,j] < avgdist && distance[o,route[i]] > avgdist
#                 println(avgdist, " ", j, " ", distance[o,j], " ", route[i], " ", distance[o,route[i]])
#                 i += 1
#             end
#             i += 1
#             if i < length(route)
#                 println(j, " ", i, " ", route[i])
#                 moveNode(route, j, indexj, i)
#             end
#         end
#     end
# end
#
# function move(route, near, index1, index2)
#     indexj = index1
#     for j in route[index1+1:index2-1]
#         indexj += 1
#         # if j in near
#             deleteat!(route, indexj)
#             println(j, " ", route[index1], " ", distance[j, route[index1]], " ", route[index2], " ", distance[j, route[index2-1]])
#             if distance[j, route[index1]] < distance[j, route[index2-1]]
#                 insert!(route, index1, j)
#                 index1 += 1
#             else
#                 insert!(route, index2, j)
#                 index2 -= 1
#             end
#             println(route)
#         # end
#     end
#     return index1, index2
# end
#
# # if slow, decide where to insert using distance from depot (avoiding zigue-zague idea)
# function moveNode(route, j, index1, index2, mind)
#     minroute = copy(originalroute)
#
#     for i in index1:index2
#         insert!(route, i, j)
#         d = getRouteDistance(route)
#         if d < mind
#             mind = d
#             minroute = copy(route)
#         end
#         deleteat!(route, i)
#     end
#     return minroute, mind
# end
#
# function moveByAngle(route, angle, index1, index2)
#     mind = getRouteDistance(route)
#     sumangle = 0
#     for j in route[2:end-1]
#         sumangle += angle[j]
#     end
#     println(sumangle/3, " moveByAngle ", 2*sumangle/3)
#
#     indexj = 0
#     for j in route[2:end-1]
#         indexj += 1
#         if angle[j] < 2*sumangle/3 && angle[j] > sumangle/3
#             deleteat!(route, indexj)
#
#             if indexj < index1
#                 route1, d = moveNode(route, j, index2, length(route), mind)
#             elseif indexj > index2
#                 route1, d = moveNode(route, j, 2, index1, mind)
#             else
#                 route1, d1 = moveNode(route, j, 2, index1, mind)
#                 route2, d2 = moveNode(route, j, index2, length(route), mind)
#                 if mind2 < d1
#                     route1 = route2
#                     d = d2
#                 else
#                     d = d1
#                 end
#             end
#
#             println(route1, " mind ", d, " ", mind)
#             if d < mind
#                 mind = d
#                 route = route1
#             end
#         end
#     end
#
#     return route
# end
#
# function splitNearFarNodes(route, avg)
#     near, far = [], []
#     for j in route[2:end-1]
#         push!(distance[o,j] <= avg ? near : far, j)
#     end
#     return near, far
# end
#
# function petalLS(originalroute, dictangle)
#     if length(originalroute) < 5
#         return originalroute
#     end
#
#     route = copy(originalroute)
#     avg = geMinMaxDistance(route, 2, length(route))
#     near, far = splitNearFarNodes(route, avg)
#     index1, index2 = findfirst(isequal(far[1]),route), findfirst(isequal(last(far)),route)
#     println(index1, " far nodes ", index2)
#
#     index1, index2 = move(route, near, index1, index2) # avoiding zig-zague in longitude for nodes between index1 and index2
#     moveUp(route, index1) # avoiding zig-zague in longitude for nodes before index1
#     moveDown(route, index2) # avoiding zig-zague in longitude for nodes after index2
#     route = getRouteDistance(route) > getRouteDistance(originalroute) ? copy(originalroute) : route
#     route = moveByAngle(route, dictangle, index1, index2)
#
#     return route # testar com poucas rotas + posprocessing retirando nós com qjrt = 0
# end
#
# function getRouteBigTruck(route)
#     while length(route) > 2 && getRouteDistance(route)/minspeed + (length(route)-1)*maxloadtime > working_hours
#         deleteat!(route, length(route)-1)
#     end
# end

# function addRoute(routes, route)
#     # println("total: ", route)
#     if length(route) > 3
#         # push!(route, o)
#         lr = length(route)
#         if lr < 11
#             aux = @elapsed route = tsp([route;o])
#             # println(route, " ", aux)
#         end
#         # cr = combinations(route[2:end], nr); # tsp([o;r;o]
#         for nr in 2:5
#             i = 2
#             while i < lr # for r in cr
#                 # nr = rand([2,2,2,3,3,3,3,3,3,4,4,4,5,5,6])
#                 index = i+1.5*nr < lr ? i+nr-1 : lr
#                 push!(routes, tsp([o;route[i:index];o])) #petalLS(route, dictangle))
#                 # println(i, " ", index, " ", routes[end])
#                 i = index+1
#             end
#         end
#     end
# end
#
# counter-clockwise, limited by travel time only
# function onePetal(routes)
#     minspeed = minimum(vspeed)[2]
#     maxloadtime = maximum(vloadtime)[2]
#     sectors, orderedDPs = getSectors()
#
#     i = 1
#     for l in sectors
#         while i < l
#             j = orderedDPs[i] #indexj
#             route = [o,j]
#             i += 1
#
#             while i < l && getRouteDistance([route;orderedDPs[i];o])/minspeed + maxloadtime + unloadtime * length(route) < working_hours
#                 push!(route, orderedDPs[i])
#                 i += 1
#             end
#             addRoute(routes, route)
#         end
#     end
# end
