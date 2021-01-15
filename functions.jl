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
    if isfile(filename)
        local aux = readdlm(filename, '$')
        for i in 1:trunc(Int, length(aux)/4)
            pie[aux[i,1],aux[i,2],aux[i,3]] = aux[i,4]
        end
        distance = floydWarshall()
        return true
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

function permute(N, l, r, all, city)
    if l == r
        push!(all, copy(N))
    else
        for i in l:r
            if l != 1 || N[i] != city
                aux = N[l]
                N[l] = N[i]
                N[i] = aux
                permute(N, l+1, r, all, city)
                aux = N[l]
                N[l] = N[i]
                N[i] = aux
            end
        end
    end
end

function tsp(N, t)
    all = Array{Array{String}}(undef,0)
    permute(N, 1, length(N)-1, all, N[length(N)-1])

    mindist = B
    best = 0
    for route in all
        dist = getRouteDistance(route, t)
        if dist < mindist
            mindist = dist
            best = route
        end
        # println(route, dist)
    end

    return best
end

function getRouteDistance(route, t)
    if length(route) > 0
        dist = 0
        u = o
        for j in route
            dist += distance[u,j,t]
            u = j
        end
        return dist
    end
    return 0
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

function updateValues(t1, t2, rent, fuel, from, to)
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
        for c in C
            for j in DP
                to[c,j,t] = from[c,j,t]
            end
        end
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
