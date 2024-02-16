include("functions.jl")
include("input.jl")
include("data.jl")
include("constructive.jl")
include("greactive.jl")
include("hc.jl")
include("petal.jl")
include("ps.jl")
include("debug.jl")

function runCplex()
    m = createModel(EXACT)
    setObjective(m)

    routesch, incumbentch, unsatisfiedpenalty, transportationcosts, eqpenalty, vset, qdict = runCH()
    for i in No
        for j in No
            for k in K
                for t in T
                    for p in P
                        set_start_value(m[:V][i,j,k,t,p], ((i,j,k,t,p) in vset) ? 1 : 0)
                    end
                end
            end
        end
    end

    elapsedtime = @elapsed optimize!(m)
    println("\n", nK, " ", total_demand, " ", objective_bound(m), " ", objective_value(m), " ", elapsedtime)
    return has_values(m)
end

function runCH()
    #global nK

    demandhc, supplyhc = copy(demand), copy(supply)
    routesch, demandhc, rent, fuel, eqpenalty, vset, qdict, loads = createSolution("onlybig", 1, nT, demandhc, supplyhc)
    transportationcosts = sum(values(rent)) + sum(values(fuel))
    unsatisfiedpenalty = sum(priorityo[c] * demandhc[c,j,t] for c in C for j in DP for t in T)
    incumbentch = unsatisfiedpenalty + transportationcosts + (EQCONSTR ? wfairness*eqpenalty : 0)

    #= reduce K set, removing unused vehicles
    unusedK = []
    for k in K
        used = false
        for t in T
            if ("theta",o,k,t,1) in vset
                used = true
                break
            end
        end
        if !used
            push!(unusedK, k)
            nK -= 1
        end
    end
    setdiff!(K, unusedK)
=#
    if DEBUG
        println(nK, " used vehicles: ", K)
        printRouting(C, No, K, T, P, vset, qdict, vload, vmodel, nothing, "CH", incumbentch, unsatisfiedpenalty)
    end

    return routesch, incumbentch, unsatisfiedpenalty, transportationcosts, eqpenalty, vset, qdict
end

function runPS(vset, qdict, incumbent, method)
    if incumbent > 0
        loads = Dict()
        pselapsedtime = @elapsed psfototal, psunsatisfied, pstransportation, pseq, inuse, unsatisfiedDemand, loads, qdict =
            proximitySearch(vset, incumbent, RERUN_LS ? B : 1)

        if DEBUG && psfototal > 1
            printRouting(C, No, K, T, P, vset, qdict, vload, vmodel, nothing, method, psfototal, psunsatisfied)
        end
        GC.gc()
        sleep(SLEEP)
        GC.gc()
        return psfototal, psunsatisfied, pstransportation, pseq
    end
    return -1*ones(4)
end

function runGRASP()
    routesch, incumbentch, vset, qdict = runCH()
    elapsedtime = @elapsed bestincumbent, bestunsatisfied, besttransportation, besteq, bestit, besttime, bestvset, bestqdict = gmain(incumbentch)
    if RUN_LS
        pselapsedtime = @elapsed psincumbent, psunsatisfied, pstransportation, pseq = runPS(bestvset, bestqdict, bestincumbent, "PS GRASP")
        if psincumbent < 1
            psincumbent, psunsatisfied, pstransportation, pseq = bestincumbent, bestunsatisfied, besttransportation, besteq
        end
        println(bestincumbent, " ", bestunsatisfied, " ", besttransportation, " ", besteq, " ",
            psincumbent, " ", psunsatisfied, " ", pstransportation, " ", pseq, " ", bestit, " ", besttime, " ", elapsedtime, " ", pselapsedtime)
    else
        println(bestincumbent, " ", bestunsatisfied, " ", bestit, " ", besttime, " ", elapsedtime)
    end
end

function main()
    chelapsedtime = @elapsed routesch, incumbentch, unsatisfiedch, transportationch, eqch, vset, qdict = runCH()
    pschelapsedtime = @elapsed psch, psunsatisfiedch, pstransportationch, pseqch =
        RUN_LS ? runPS(vset, qdict, incumbentch, "PS CH") : (incumbentch, unsatisfiedch, transportationch, eqch)
    if psch < 1
        psch, psunsatisfiedch, pstransportationch, pseqch = incumbentch, unsatisfiedch, transportationch, eqch
    end

    routespetal = Dict()
    #=petalonlyelapsedtime = @elapsed incumbentpetal, unsatisfiedpetalonly, transportationpetalonly, eqpetalonly,
        vset, qdict, loads = runPetal(routespetal, "Petal only")
    pspetalonlyelapsedtime = @elapsed pspetalonly, psunsatisfiedpetalonly, pstransportationpetalonly, pseqpetalonly =
        false ? runPS(vset, qdict, incumbentpetal, "PS Petal only") : (incumbentpetal, unsatisfiedpetalonly, transportationpetalonly, eqpetalonly)
    if pspetalonly < 1
        pspetalonly, psunsatisfiedpetalonly, pstransportationpetalonly, pseqpetalonly = incumbentpetal, unsatisfiedpetalonly, transportationpetalonly, eqpetalonly
    end=#

    if FIRSTCHOICE == BYPRIORITY
    petalelapsedtime = @elapsed incumbent, unsatisfiedpetal, transportationpetal, eqpetal,
        vset, qdict, loads = runPetal(routespetal, "Petal", routesch)
    pspetalelapsedtime = @elapsed pspetal, psunsatisfiedpetal, pstransportationpetal, pseqpetal =
        RUN_LS ? runPS(vset, qdict, incumbent, "PS Petal") : (incumbent, unsatisfiedpetal, transportationpetal, eqpetal)
    if pspetal < 1
        pspetal, psunsatisfiedpetal, pstransportationpetal, pseqpetal = incumbent, unsatisfiedpetal, transportationpetal, eqpetal
    end
    end

    #petalelapsedtime1 = @elapsed incumbent1, unsatisfiedpetal1, transportationpetal1, eqpetal1,
    #    vset, qdict, loads = runPetal(routespetal, "Petal", routesch, 2)

    #println("\n", instance, " ", nK, " ", incumbentch, " ", unsatisfiedch, " ", transportationch, " ", eqch)
    print("\n", instance, " ", nK, " ", #sum(length(r) for r in values(routesch)), " ", minimum(length(routespetal[t]) for t in T), " ", maximum(length(routespetal[t]) for t in T, " "),
        incumbentch, " ", unsatisfiedch, " ", transportationch, " ", eqch, " ")
	#incumbentpetal, " ", unsatisfiedpetalonly, " ", transportationpetalonly, " ", eqpetalonly, " ",
    if FIRSTCHOICE == BYPRIORITY
        print(incumbent, " ", unsatisfiedpetal, " ", transportationpetal, " ", eqpetal, " ",
        #incumbent1, " ", unsatisfiedpetal1, " ", transportationpetal1, " ", eqpetal1, " ",
        #psch, " ", psunsatisfiedch, " ", pstransportationch, " ", pseqch, " ",
        #pspetalonly, " ", psunsatisfiedpetalonly, " ", pstransportationpetalonly, " ", pseqpetalonly, " ",
        #pspetal, " ", psunsatisfiedpetal, " ", pstransportationpetal, " ", pseqpetal, " ",
        chelapsedtime, " ", petalelapsedtime) end
    println()
end

if METHOD == EXACT
    runCplex()
elseif METHOD == GRASP
    runGRASP()
else
    main()
end
