include("functions.jl")
include("input.jl")
include("data.jl")
include("constructive.jl")
include("ps.jl")
include("psbound.jl")
include("lineq.jl")
include("debug.jl")

hcrent, hcfuel, hcfo, hceqsum, hceqmean, hceqmax, hcelapsedtime = Array{Float64}(undef, 6), Array{Float64}(undef, 6), Array{Float64}(undef, 6), Array{Float64}(undef, 6), Array{Float64}(undef, 6), Array{Float64}(undef, 6), Array{Float64}(undef, 6)
psrent, psfuel, psfo, pseqsum, pseqmean, pseqmax, pselapsedtime = Array{Float64}(undef, 6), Array{Float64}(undef, 6), Array{Float64}(undef, 6), Array{Float64}(undef, 6), Array{Float64}(undef, 6), Array{Float64}(undef, 6), Array{Float64}(undef, 6)
linfo, lineqsum, withinBound = Array{Float64}(undef, 6), Array{Float64}(undef, 6), Array{Bool}(undef, 6)
nK2 = Array{Integer}(undef, 6)

if METHOD == EXACT
    tstart = time()
    vset, qdict = Set(), Dict()
    unsatisfiedDemand = Dict()
    for c in C
        for j in DP
            for t in T
                unsatisfiedDemand[c,j,t] = 0
            end
        end
    end
    m = createModel(C, DP, N, K, 1, nT, supply, demand, TIME_LIMIT)
    #=fo, rent, fuel, eq, eqmax, gap = solveModel(m, C, DP, N, K, T, 0, vset, qdict, unsatisfiedDemand)
    elapsedtime = time()-tstart
    if rent > 0
        printRouting(C, No, K, T, P, vset, qdict, vload, vmodel, unsatisfiedDemand, TIME_LIMIT)
    end
    println(nK, " ", total_demand, " ", fo, " ", sum(priority[c] * unsatisfiedDemand[c,j,t] for c in C for j in DP for t in Tend[c]), " ", rent, " ", fuel, " ", eq, " ", eqmax, " ", elapsedtime, " ", termination_status(m))=#
else
    if KNOWFUTURE
        for t in 2:nT-2
            getDataUpdate(t, 0)
        end
    end

    methods = ["onlybig"] # "notrip", "1:nK", "nK:1", "onlybig", "onlysmall", "both"
    besthc, bestps = "", ""
    fohc, fops, minK = B, B, B

    for im in 1:length(methods)
        demandhc = copy(demand)
        supplyhc = copy(supply)
        demandhc, loads, rent, fuel, hceqsum[im], hceqmean[im], hceqmax[im], vset, qdict, hcelapsedtime[im] = createSolution(methods[im], 1, nT, demandhc, supplyhc)
        hcfo[im] = sum(priority[c] * demandhc[c,j,t] for c in C for j in DP for t in T)
        hcrent[im] = sum(values(rent))
        hcfuel[im] = sum(values(fuel))
        aux = hcfo[im] + hcrent[im] + hcfuel[im] + wfairness*hceqsum[im]
        if aux < fohc
            global besthc = methods[im]
            global fohc = aux
        end

        if DEBUG
            printRouting(C, No, K, T, P, vset, qdict, vload, vmodel, demandhc, hcelapsedtime[im])
        end

        if METHOD == PS
            supplyps = copy(supply)
            demandps, rentps, fuelps, pseqsum[im], pseqmean[im], pseqmax[im], pselapsedtime[im], vset, qdict, loads, nK2[im] = proximitySearch(3, nT, supplyps, demandhc, loads, rent, fuel, vset, qdict)
            GC.gc()
    	    sleep(SLEEP)
            linfo[im], lineqsum[im], qdictlin = linearEq(loads, vset, qdict, psfo[im], 0, 0)

            psfo[im] = sum(priority[c] * demandps[c,j,t] for c in C for j in DP for t in T)
            psrent[im] = sum(values(rentps))
            psfuel[im] = sum(values(fuelps))
            aux = psfo[im] + psrent[im] + psfuel[im] + wfairness*pseqsum[im]
            if aux < fops
                global bestps = methods[im]
                global fops = aux
            end

            global minK = min(minK, nK2[im])

            # withinBound[im] = proximitySearchB(vset, aux)
        else
            linfo[im], lineqsum[im], qdictlin = linearEq(loads, vset, qdict, hcfo[im], 0, 0)
        end

        GC.gc()
	    sleep(SLEEP)
    end

    for im in 1:length(methods)
        print("\n", methods[im], " ", hcrent[im], " ", hcfuel[im], " ", hcfo[im], " ", hceqsum[im], " ", hceqmean[im], " ", hceqmax[im], " ", hcelapsedtime[im])
        if METHOD == PS print(" ", psrent[im], " ", psfuel[im], " ", psfo[im], " ", pseqsum[im], " ", pseqmean[im], " ", pseqmax[im], " ", pselapsedtime[im], " ", withinBound[im]) end
        print(" ", linfo[im], " ", lineqsum[im])
    end
    println("\n", nK, " ", minK, " ", total_demand, " ", fohc, " ", fops, " ", besthc, " ", bestps)
end
