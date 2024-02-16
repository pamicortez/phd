include("functions.jl")
include("input.jl")
include("data.jl")
include("constructive.jl")
include("hc.jl")
include("ps.jl")
include("debug.jl")

function print_debug(supplyh, loadshc, loadsps, demandhc, unsatdemand)
    for c in C
        t1 = 1
        for t in Tend[c]
            println("\"", c, "\" ", t, " ", sum(supply[c,t2] for t2 in t1:t), " ", sum(supplyh[c,t2] for t2 in t1:t))
            for j in DP
                println("\"", c, "\" \"", j, "\" ", t, " ", sum(demand[c,j,t2] for t2 in t1:t), " ",
                    sum(loadshc[c,j,t2] for t2 in t1:t), " ", sum(loadsps[c,j,t2] for t2 in t1:t), " ",
                    sum(demandhc[c,j,t2] for t2 in t1:t), " ", unsatdemand[c,j,t])
            end
            t1 = t+1
        end
    end
end

# deterministic evaluation
function deterministic(loadshc, loadsps)
    for t in T
        getDataUpdate(t, Dict())
    end
    demandhc = copy(demand)
    unsatdemand = copy(demand)

    for t in T
        updateResidualDemand(C, DP, t, 1, demandhc, loadshc)
        if loadsps != -1 updateResidualDemand(C, DP, t, 1, unsatdemand, loadsps) end
    end

    dhc, eqhc = sum(priority[c] * demandhc[c,j,t] for c in C for j in DP for t in T), sum(values(getEquity(demandhc)))
    dps = loadsps != -1 ? sum(priority[c] * unsatdemand[c,j,t] for c in C for j in DP for t in T) : dhc
    eqps = loadsps != -1 ? sum(values(getEquity(unsatdemand))) : eqhc

    return dhc, eqhc, dps, eqps
end

function antecipative(method)
    if KNOWFUTURE
        for t in T
            getDataUpdate(t, Dict())
        end
    end
    supplyhc = copy(supply)
    demandhc = copy(demand)

    routes, demandhc, rent, fuel, hceqsum, vset, qdict, loadshc = createSolution(method, 1, nT, demandhc, supplyhc)
    hcfo = sum(priority[c] * demandhc[c,j,t] for c in C for j in DP for t in T)
    r, f = getValues(1, nT, rent, fuel)
    hcrent = r
    hcfuel = f
    hcfototal = hcrent + hcfuel + hcfo + (EQCONSTR ? wfairness*hceqsum : 0)
    if DEBUG
        printRouting(C, No, K, T, P, vset, qdict, vload, vmodel, demandhc, "CH", hcfo, hcfototal)
    end

    pselapsedtime = @elapsed fo, psfo, psfuel, pseqsum, inuse,
        unsatdemand, loadsps, qdict = proximitySearch(vset, hcfototal, 1)


    if fo < 1
    #     if DEBUG
    #         printRouting(C, No, K, T, P, vset, qdict, vload, vmodel, unsatdemand, "PS", pselapsedtime, fo)
    #         print_debug(supplyhc, loadshc, loadsps, demandhc, unsatdemand)
    #     end
    # else
        psfo, psrent, psfuel, pseqsum = hcfo, hcrent, hcfuel, hceqsum
        fo = psfo + psrent + psfuel + (EQCONSTR ? wfairness*pseqsum : 0)
    end

    if !KNOWFUTURE
        hcfo, hceqsum, psfo, pseqsum = deterministic(loadshc, loadsps)
    end

    print("\n", method, " ", hcfo, " ", hcrent + hcfuel, " ", hceqsum)
    println(" ", psfo, " ", psrent + psfuel, " ", pseqsum, " ", pselapsedtime)

    return hcfototal, fo, pselapsedtime
end

function dynamic(method)
    hcfototal, hcfo, hcrent, hcfuel, rent, fuel, hceqsum, hceqmean, hceqmax, pselapsedtime,
        fo, psrent, psfuel, psfo, pseqsum = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    v, q, vset, qdict, loadshc, loadsps, unsatdemand = nothing, nothing, nothing, nothing, nothing, nothing, nothing
    supplyhc = copy(supply)
    demandhc = copy(demand)

    for t in T
        if t > 1
            updated = getDataUpdate(t, q)
            supplyhc, demandhc = getResidual(t, loadsps)
        end

        if t == 1 || updated
            routes, demandhc, rent, fuel, hceqsum, vset, qdict, loadshc = createSolution(method, t, nT, demandhc, supplyhc)
            hcfo = sum(priority[c] * demandhc[c,j,t2] for c in C for j in DP for t2 in T)
            hceqsum = sum(values(getEquity(demandhc)))#, Statistics.mean(values(eq)), maximum(values(eq))

            if DEBUG
                printRouting(C, No, K, T, P, vset, qdict, vload, vmodel, demandhc, method, hcfo, nothing, t)
            end
        end
        r, f = getValues(t, t, rent, fuel)
        hcrent += r
        hcfuel += f
        hcfototal = hcrent + hcfuel + sum(getValues(t+1, nT, rent, fuel)) + hcfo + (EQCONSTR ? wfairness*hceqsum : 0)

        if t == 1
            v = copy(vset)
            q = copy(qdict)
        else
            copySolution(t:nT, vset, qdict, v, q)
        end

        m = createModel()
        pselapsedtime += @elapsed incumbent, pseqsum2, unsatdemand2, loadsps2, q2, psrent2, psfuel2 =
            proximitySearch(m, t, v, hcfototal, pselapsedtime, q)
        if incumbent > 0
            pseqsum, unsatdemand, loadsps, q, psrent, psfuel = pseqsum2, unsatdemand2, loadsps2, q2, psrent2, psfuel2
            psfo = sum(priority[c] * unsatdemand[c,j,t2] for c in C for j in DP for t2 in Tend[c])
            fo = psfo + psrent + psfuel + (EQCONSTR ? wfairness*pseqsum : 0)
            psfuel += psrent
        else
            fo, psfo, psfuel, pseqsum, inuse, unsatdemand, loadsps, q = solveFixedModel(m, v, q, t)
        end
        m = nothing
        if DEBUG
            printRouting(C, No, K, T, P, v, q, vload, vmodel, unsatdemand, "PS "*method, fo, pselapsedtime, t)
            print_debug(supplyhc, loadshc, loadsps, demandhc, unsatdemand)
        end
    end

    # print_debug(loadshc, loadsps, demandhc, unsatdemand)
    print("\n", method, " ", hcfo, " ", hcrent + hcfuel, " ", hceqsum)
    println(" ", psfo, " ", psfuel, " ", pseqsum, " ", pselapsedtime)

    return hcfototal, fo, pselapsedtime
end

function main()
    method = "notrip"
    # methods = ["notrip", "1:nK", "nK:1", "onlybig", "onlysmall", "both"]
    # besthc, bestps = "", ""
    # fohc, fops = B, B

    # for method in methods
        hcfototal, psfototal, elapsedtime = DYNAMIC_ANALYSIS ? dynamic(method) : antecipative(method)

        # if hcfototal > 0 && fohc > hcfototal
        #     fohc = hcfototal
        #     besthc = method
        # end
        # if psfototal > 0 && fops > psfototal
        #     fops = psfototal
        #     bestps = method
        # end

        # global supply, demand = getSupplyDemand(psupply)
        # GC.gc()
        # sleep(SLEEP)
    # end
    # println(nK, " ", total_demand, " ", fohc, " ", fops, " ", besthc, " ", bestps, " ", elapsedtime/6)
    println(nK, " ", total_demand, " ", hcfototal, " ", psfototal, " ", elapsedtime)
end

main()
