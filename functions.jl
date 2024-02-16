using Printf

# para criar matrix para cada t, multiplicando distancia por pi_ijt
function floydWarshall()
    aux = Array{String,1}(undef, n)
    dist0 = Array{Float64,2}(undef, n, n)
    dist = Array{Float64,3}(undef, n, n, nT)
    distance = Dict{Any,Float64}()

    ii = jj = 1
    for i in 1:645
        if cities[i] in N
            for j in 1:645
                if cities[j] in N
                    dist0[ii,jj] = auxDistance[i,j]
                    jj += 1
                end
            end
            aux[ii] = cities[i]
            ii += 1
            jj = 1
        end
    end

    for i in 1:n
        for j in 1:n
            for t in T
                dist[i,j,t] = i == j ? 0 : (dist0[i,j] * (pie[aux[i],aux[j],t] > 0.9 ? pie[aux[i],aux[j],t] : B))
            end
        end
    end

    for t in T
        for k in 1:n
            for i in 1:n
                for j in 1:n
                    if dist[i,j,t] > dist[i,k,t] + dist[k,j,t]
                        dist[i,j,t] = dist[i,k,t] + dist[k,j,t]
                    end
                end
            end
        end
    end

    for t in T
        for i in 1:n
            for j in 1:n
                distance[aux[i],aux[j],t] = dist[i,j,t]
            end
        end
    end

    return distance #maximum(dist),
end

function readPie(filename)
    changed = false
    if isfile(filename)
        local aux = readdlm(filename, '$')
        for i in 1:trunc(Int, length(aux)/4)
            if !((aux[i,1],aux[i,2],aux[i,3]) in keys(pie)) || pie[aux[i,1],aux[i,2],aux[i,3]] != aux[i,4]
                pie[aux[i,1],aux[i,2],aux[i,3]] = aux[i,4]
                changed = true
            end
        end
        if changed
            global distance = floydWarshall()
        end
        return changed
    end
    return false
end

function minGreaterZero(d::Matrix)
    minval = typemax(Float64)

    for i in d
        if i > 0 && i < minval
            minval = i
        end
    end

    minval
end

function minDistance(d::Dict{Any,Float64})
    minval = typemax(Float64)

    for t in T
        for j in DP
            if d[o,j,t] < minval
                minval = d[o,j,t]
            end
        end
    end

    return minval
end

function nearestNeighbours(t, u, prio)
    i = 0
    aux = Array{Float64}(undef, nDP)
    for j in DP
        i += 1
        # prioj = sum(prio[c,j] for c in C)
        aux[i] = distance[u,j,t]
    end
    return DP[sort([1:nDP;], by=i->(aux[i]))]
end

function halfpermute(original, result)
    for i in 2:length(original)-1
        for j in i+1:length(original)-1
            N = copy(original)
            aux = N[2]; N[2] = N[i]; N[i] = aux
            aux = N[end-1]; N[end-1] = N[j]; N[j] = aux
            permute(N, 3, length(N)-2, result)
        end
    end
end

function permute(N, l, r, result)
    if l == r
        push!(result, copy(N))
    else
        for i in l:r
            aux = N[l]
            N[l] = N[i]
            N[i] = aux
            permute(N, l+1, r, result)
            aux = N[l]
            N[l] = N[i]
            N[i] = aux
        end
    end
end

function tsp(N, t)
    if length(N) <= 4
        return N
    end

    result = Array{Array{String}}(undef,0)
    halfpermute(N, result)
    # println("tsp: ", length(result))

    mindist = B
    best = 0
    for route in result
        dist = getRouteDistance(route, t)
        if dist < mindist
            mindist = dist
            best = route
        end
        # println(length(result), route, dist)
    end

    return best
end

function getRouteDistance(route, t)
    dist = 0
    if length(route) > 0
        u = o
        for j in route
            dist += distance[u,j,t]
            u = j
        end
    end
    return dist
end

function getRouteDistance(route)
    dist = 0
    if length(route) > 0
        u = route[1]
        for j in route[2:end]
            dist += distance[u,j]
            u = j
        end
    end
    return dist
end

function getOrderByNearest(route, t)
    result = [o]
    u = o
    insolution = [false for i in 1:length(route)]
    for i in 2:length(route)-1
        nearest = B
        nearestj = 0
        for j in 2:length(route)-1
            if !insolution[j] && distance[u,route[j],t] < nearest
                nearest = distance[u,route[j],t]
                nearestj = j
            end
        end
        push!(result, route[nearestj])
        insolution[nearestj] = true
        u = route[nearestj]
        # println(j, " ", nearestj, " ", insolution)
    end
    push!(result, o)
    # println(route)
    # println(result)
    return result
end

function getMinTripTime(k, t)
    minTripTime = B
    for j in DP
        aux = (distance[o,j,t] + distance[j,o,t])/vspeed[vmodel[k]] + vloadtime[vmodel[k]]
        if minTripTime > aux
            minTripTime = aux
        end
    end
    return minTripTime
end

# roads condition - na verdade as estradas bloqueadas sao as mais faceis de resolver: basta retirar arvores. 1.2 e 1.5 sao danos estruturais
function setPi(pie, i, j, value)
    for t in T
        r = rand()
        if (value == 1.5 && r/t <= 0.047) || (value == 1.2 && r/t <= 0.1) || (value == 0 && r/t <= 0.167)
            value = 1
        end

        pie[i,j,t] = value
        pie[j,i,t] = value
    end
end

function getPi()
    # initialisation: roads condition
    pie = Dict{Any,Float32}()
    for i in N
        for j in N
            for t in T
                pie[i,j,t] = 1
                pie[j,i,t] = 1
            end
        end
    end

    epicentres = unique(rand(1:nDP, 3))
    for e in epicentres
        epicentre = DP[e]
        for i in N
            if distance[epicentre, i] > ZERO
                r = rand()
                if r < 0.7
                    setPi(pie, epicentre, i, r < 0.05 ? 0 : (r < 0.3 ? 1.5 : 1.2))
                end
            end
        end
    end

    return pie
end

function getSupplyDemand()
    supply = Dict{Any,Float64}()
    demand = Dict{Any,Float64}()
    for t in T
        supply["food kit",t]     = (t-1) % 5 == 0 ? trunc(Integer, total_demand) : 0
        supply["cleaning kit",t] = t == 1 || t == 11 ? trunc(Integer, total_demand/4) : 0
        supply["mattress",t]     = t == 1 ? trunc(Integer, total_demand) : 0

        for (j,pop) in peopleDP
            demand["food kit",j,t]     = (t-1) % 5 == 0 ? pop : 0
            demand["cleaning kit",j,t] = t == 1 || t == 11 ? trunc(Integer, pop/4) : 0
            demand["mattress",j,t]     = t == 1 ? pop : 0
        end
    end
    return supply, demand
end

function getValues(t1, t2, rent, fuel)
    r = f = 0
    for t in t1:t2
        if isa(findmin(rent)[2], Number)
            r += rent[t]
        else
            for k in K
                r += rent[k,t]
            end
        end
        f += fuel[t]
    end
    return r, f
end

function copySolution(C, No, K, T, fromv, fromq, tov, toq)
    for k in K
        for t in T
            for p in P
                for i in No
                    for j in No
                        if (i,j,k,t,p) in fromv
                            push!(tov, (i,j,k,t,p))
                            # print("V ", k, " ", t, " ", p, " (", i, " -> ", j, ")")
                            if i != "theta" && j != "theta" && j != o
                                for c in C
                                    toq[c,i,j,k,t,p] = fromq[c,i,j,k,t,p]
                                    # @printf("\t%s = %.2f\t ", c, toq[c,i2,i,j,k,t,p])
                                end
                            end
                            # println()
                        end
                    end
                end
            end
        end
    end
end

function copySolution(T, fromv, fromq, tov, toq)
    for k in K
        for t in T
            for p in P
                for i in No
                    for j in No
                        if (i,j,k,t,p) in tov
                            delete!(tov, (i,j,k,t,p))
                        end
                        if (i,j,k,t,p) in fromv
                            push!(tov, (i,j,k,t,p))
                        end
                        if i != "theta" && j != "theta" && j != o
                            for c in C
                                if (c,i,j,k,t,p) in keys(toq)
                                    toq[c,i,j,k,t,p] = 0
                                end
                                if (c,i,j,k,t,p) in keys(fromq)
                                    toq[c,i,j,k,t,p] = fromq[c,i,j,k,t,p]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    # for value in values(vset)
    #     if value[4] == t
    #         push!(v, value)
    #     end
    # end
    # for key in keys(qdict)
    #     if key[5] == t
    #         q[key] = qdict[key]
    #     end
    # end
end

function get_available_time(k, t, vset, p2=0)
    if !vehicle[k,t]
        return 0
    end

    auxdistance = auxp = 0
    aux = collect(2:nP)
    localP = deleteat!(aux, findall(x->x==p2,aux))

    for p in localP
        isp_used = false
        for i in N
            for j in filter(e -> e != i, N)
                if isVar(i,j,p) && (i,j,k,t,p) in vset
                    isp_used = true
                    auxdistance += distance[i,j,t]
                end
            end
        end
        if isp_used
            auxp += 1
        end
    end

    return working_hours - auxdistance/vspeed[vmodel[k]] - auxp*vloadtime[vmodel[k]]
end

function get_available_load(k,t,p, originalroute, route, qdict)
    unchangedDPs = filter(e->e!=o, intersect(originalroute,route))
    return vload[vmodel[k]] - (length(unchangedDPs) > 0 ? sum(weight[c] *
        (qdict[c,originalroute[findall(x->x==j, originalroute)[1]-1],j,k,t,p] -
        get(qdict,(c,j,originalroute[findall(x->x==j, originalroute)[1]+1],k,t,p), 0)) for c in C for j in unchangedDPs) : 0)
end

function get_route(k,t,p, vset)
    route = [o]
    arcs = []
    nodes = copy(N)

    u = o
    j = nothing
    while j != o
        for i in 1:length(nodes)
            j = nodes[i]
            if (u,j,k,t,p) in vset && !((u,j) in arcs)
                push!(route, j)
                push!(arcs, (u,j))
                u = j
                break
            end
        end
    end

    # println("get_route ", k, " ", t, " ", p, " ", route)
    return route
end

function get_route_load(k,t,p, route, qdict)
    routeload = Dict()
    # routetime = vloadtime[vmodel[k]]
    if length(route) > 2
        route = copy(route)
        u = popfirst!(route)
        i = popfirst!(route)
        routeload[u] = 0
        while length(route) > 0
            j = popfirst!(route)
            routeload[i] = sum(priority[c]*(qdict[c,u,i,k,t,p] - get(qdict, (c,i,j,k,t,p), 0)) for c in C)
            # routetime += distance[u,i,t]
            u = i
            i = j
        end
        # routetime += distance[u,i,t]
        # routetime /= vspeed[vmodel[k]]
    end
    return routeload #, routetime
end

function get_routes_loads(k, t, vset, qdict)
    routes, routeloads = [], []

    for p in 2:nP
        route = get_route(k,t,p, vset)
        if length(route) < 2
            break
        end
        push!(routes, route)
        push!(routeloads, sum(values(get_route_load(k,t,p, route, qdict))))
    end

    return routes, routeloads
end
