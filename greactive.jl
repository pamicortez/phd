function addNodeG(k, t, u, prio, dday, availableTime, sortedDPprio, visited)
    return u == o ? getDP(o, k, t, availableTime, visited, dday, rand(1:2) == 1 ? prio : dday) :
            getDP(u, k, t, availableTime, visited, dday, 0)
end

function getDP(u, k, t, availableTime, visited, dday, matrix)
    g = Dict()
    for i in DP
        if !(i in visited) && (distance[u,i,t] + distance[i,o,t]) <= vspeed[vmodel[k]]*availableTime &&
            sum(dday[c,i] for c in C) > ZERO
            g[i] = length(matrix) > 1 ? -sum(matrix[c,i] for c in C) : distance[u,i,t]
        end
    end
    # println(matrix, "\n\n", g)

    j = getFromRCL(g)
    # println(j)
    return j == -1 ? "" : j
end

function getBestTruckG(t, originalk, route, routetime, routeweight, availableTime, rent)
    bestvmodel = getBestVmodel(t, t, originalk, route, routeweight, availableTime)
    if bestvmodel == ""
        bestvmodel = vmodel[originalk]
    end

    g = Dict()
    for k2 in nK:-1:1
        k = K[k2]
        if vmodel[k] == bestvmodel && routetime < availableTime[k,t]
            g[k] = (rent[k,t] < 1 ? vcost[vmodel[k]] : 0) + availableTime[k,t]
        end
    end
    k = getFromRCL(g)

    # println("k ", length(candidates))
    return k == -1 ? originalk : k
end

function getBestPeriod(tm, supplyhc, demandhc, vset, qdict, rent, fuel, loads, lastp, availableTime)
    g, aux = Dict(), Dict()
    for t in tm:tend[C[1],tm]
        fo, k, route = createTentativeRoute(t, supplyhc, demandhc, vset, qdict, rent, fuel, loads, availableTime, lastp)
        if k != "" && length(route) > 1
            g[t] = fo
            aux[t] = k
        end
    end

    t = getFromRCL(g)
    if t != -1
        return t, aux[t]
    end
    return -1, -1
end

# RCL = {i ∈ C | g(i) ≤ gmin + α(gmax − gmin)};
# α = 0 corresponds to a pure greedy algorithm, while α = 1 is equivalent to a random construction
function getFromRCL(g)
    candidates = []
    if length(g) > 0
        gmin = minimum(values(g))
        gmax = maximum(values(g))

        for i in keys(g)
            # println(g[i], " ", gmin + alpha * (gmax - gmin))
            if g[i] <= gmin + alpha * (gmax - gmin) + ZERO
                push!(candidates, i)
            end
        end
    end
    # println(candidates)
    return length(candidates) == 0 ? -1 : candidates[rand(1:length(candidates))]
end

# let Ai be the average value of all solutions found using α = αi
# in my case, I calculate this average when needed by doing: (A[i]/nGrasp[i])
function updateA(fo, A, nGrasp)
    i = Integer(alpha*10)
    A[i] += fo
    nGrasp[i] += 1
end

function updateProb(A, nGrasp, probG)
    q = Array{Float64}(undef, 10)
    for i in 1:10
        q[i] = A[i] != 0 ? nGrasp[i]/A[i] : 0 # incumbent / (A[i]/nGrasp[i]) # why qi = incumbent/Ai ? just so q[i] is inversely proportional to A[i]?
    end
    totalq = sum(q)
    for i in 1:10
        probG[i] = q[i]/totalq
    end
    if DEBUG println(probG) end
end

function getAlpha(it, A, nGrasp, probG, alphaList)
    if it%ITUPDATE == 0
        updateProb(A, nGrasp, probG)
    end

    r = rand()
    i = cumulative = 0
    while cumulative < r
        i += 1
        cumulative += probG[i]
    end

    return alphaList[i]
end

function gmain(incumbentch)
    start = time()
    bestincumbent = incumbentch
    bestincumbent2 = B
    bestit = besttime = bestunsatisfied = besttransportation = besteq = -1
    bestvset, bestqdict = Set(), Dict()

    it = elapsedtime = 0
    alphaList = collect(0.1:0.1:1)
    A = zeros(10) # solution quality per alpha
    nGrasp = zeros(Int, 10) # number of solutions per alpha
    probG = fill(0.1, 10) # uniform

    while elapsedtime < 1800
        it += 1
        global alpha = getAlpha(it, A, nGrasp, probG, alphaList)
        demandhc = copy(demand)
        supplyhc = copy(supply)

        vset, qdict = Set(), Dict()
        loads = Dict()
        rent, fuel = Dict{Any, Float64}(), Dict{Integer, Float64}()
        availableTime, lastp = Dict{Any, Float64}(), Dict{Any, Integer}()
        setVars(loads, fuel, rent, availableTime, lastp)

        t = 1
        while t <= nT
            while true
                bestt, bestk = getBestPeriod(t, supplyhc, demandhc, vset, qdict, rent, fuel, loads, lastp, availableTime)
                if bestt == -1
                    break
                end
                saveIncumbent(bestk, bestt, supplyhc, demandhc, vset, qdict, rent, fuel, loads, lastp, availableTime)
            end
            t = tend[C[1],t]+1
        end

        unsatisfied = sum(priority[c] * demandhc[c,j,t] for c in C for j in DP for t in T)
        rentg = sum(values(rent))
        fuelg = sum(values(fuel))
        eq = sum(values(getEquity(demandhc)))
        incumbent = unsatisfied + rentg + fuelg + (EQCONSTR ? wfairness * eq : 0)

        # print(it, " ", round(alpha; digits = 2), " ", round(incumbent; digits = 1), " ", unsatisfied, " ")#, " ", round(rentg; digits = 1), " ", round(fuelg; digits = 1))
        updateA(incumbent, A, nGrasp)

        elapsedtime = time() - start
        if incumbent+1 < bestincumbent
            bestincumbent = incumbent
            bestunsatisfied = unsatisfied
            besttransportation = rentg + fuelg
            besteq = eq

            bestit = it
            besttime = elapsedtime
            bestvset = copy(vset)
            bestqdict = copy(qdict)
            # psch = runPS(vset, qdict, incumbent, "PS GRASP")
        end

        if DEBUG && it%250 == 0
            if bestincumbent2 > bestincumbent
                bestincumbent2 = bestincumbent
            end
            println(it, " ", round(bestincumbent2; digits = 1), " ", bestunsatisfied, " ", round(elapsedtime; digits = 2))
        end
    end

    return bestincumbent, bestunsatisfied, besttransportation, besteq, bestit, besttime, bestvset, bestqdict
end
