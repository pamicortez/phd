function createModel(C, DP, N, K, t1, t2, supply, demand, timelimit)
    Tp = T[t1:t2]
    # println(C, K, Tp, timelimit)

    mipdisplay = 0
    emphasis = 1 # CPX_MIPEMPHASIS_FEASIBILITY
    # JuMP.set_silent(m)
    if METHOD == EXACT
        timelimit = 86400
        mipdisplay = 2
        emphasis = 3 # CPX_MIPEMPHASIS_BESTBOUND
    end
    m = Model(CPLEX.Optimizer)
    set_optimizer_attribute(m, "CPX_PARAM_TILIM", timelimit)
    #=set_optimizer_attribute(m, "CPX_PARAM_FLOWCOVERS", 1)
    set_optimizer_attribute(m, "CPX_PARAM_MIRCUTS", 1)
    set_optimizer_attribute(m, "CPX_PARAM_BTTOL", 0.1)
    set_optimizer_attribute(m, "CPX_PARAM_HEURFREQ", 100)
    set_optimizer_attribute(m, "CPX_PARAM_PROBE", -1)=#
    set_optimizer_attribute(m, "CPXPARAM_Emphasis_MIP", emphasis)
    set_optimizer_attribute(m, "CPX_PARAM_THREADS", NTHREADS)
    set_optimizer_attribute(m, "CPX_PARAM_EPGAP",0.01)
    set_optimizer_attribute(m, "CPX_PARAM_MIPDISPLAY", mipdisplay)
    set_optimizer_attribute(m, "CPX_PARAM_NUMERICALEMPHASIS", 1)
    set_optimizer_attribute(m, "CPX_PARAM_EPINT", 0.0000001)

    @variable(m, V[["theta";N],["theta";N],K,Tp,P], Bin)
    @variable(m, Q[C,N,N,K,Tp,P] >= 0)
    @variable(m, deltaQ[C,DP,Tp] >= 0)

    if METHOD == PS
        Texp = []
        for c in C
            Texp = union(Texp,[tend[c,t1],tend[c,t2]])
        end
        @variable(m, dexp[C,DP,Texp] >= 0)

        @constraint(m, demandexpire[c in C, j in DP, t in unique([tend[c,t1],tend[c,t2]])],
            sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:min(t,t2)) -
            sum(deltaQ[c,j,tau] for tau in max(div(t-1,threshold[c]+1)*(threshold[c]+1)+1,t1):min(t,t2)) <= dexp[c,j,t])
        @constraint(m, nonantecipativity[c in C, j in DP, t in setdiff(Tp, unique([tend[c,t1],tend[c,t2]]))],
            sum(deltaQ[c,j,tau] for tau in max(div(t-1,threshold[c]+1)*(threshold[c]+1)+1,t1):t)
            <= sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t))

        if PSEQ
            @variable(m, 0 <= er[C,Texp] <= 1)
            @objective(m, Min, sum(priority[c] * dexp[c,j,t] for c in C, j in DP, t in Texp)
                + wfairness * sum(er[c,t] for c in C, t in Texp)
                + sum(vcost[vmodel[k]] * V["theta",o,k,t,1] for k in K, t in Tp)
                + sum(vfuel[vmodel[k]] * distance[i,j,t] * V[i,j,k,t,p] for i in N, j in N, k in K, t in Tp, p in P))

            for c in C
                for j in DP
                    for j2 in DP
                        for t in unique([tend[c,t1],tend[c,t2]])
                            if sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:min(t,t2)) > ZERO &&
                                    sum(demand[c,j2,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:min(t,t2)) > ZERO
                                @constraint(m, dexp[c,j,t]/sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:min(t,t2)) -
                                    dexp[c,j2,t]/sum(demand[c,j2,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:min(t,t2)) <= er[c,t])
                            end
                        end
                    end
                end
            end
       else
            @objective(m, Min, sum(priority[c] * dexp[c,j,t] for c in C, j in DP, t in Texp)
                + sum(vcost[vmodel[k]] * V["theta",o,k,t,1] for k in K, t in Tp)
                + sum(vfuel[vmodel[k]] * distance[i,j,t] * V[i,j,k,t,p] for i in N, j in N, k in K, t in Tp, p in P))
        end

        # for c in C
        #     for j in DP
        #         for t in unique([tend[c,t1],tend[c,t2]])
        #             println(c, " ", j, " ", t, " ", sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:min(t,t2)))
        #         end
        #     end
        # end
        # for c in C
        #     for j in DP
        #         for t in setdiff(Tp, unique([tend[c,t1],tend[c,t2]]))
        #             println(c, " ", j, " ", t, " ", sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t))
        #         end
        #     end
        # end
    else
        @variable(m, dexp[C,DP,Tendall] >= 0)
        @variable(m, 0 <= er[C,Tendall] <= 1)

        @objective(m, Min, sum(priority[c] * dexp[c,j,t] for c in C, j in DP, t in Tend[c])
            + wfairness * sum(er[c,t] for c in C, t in Tend[c])
            + sum(vcost[vmodel[k]] * V["theta",o,k,t,1] for k in K, t in T)
            + sum(vfuel[vmodel[k]] * distance[i,j,t] * V[i,j,k,t,p] for i in N, j in N, k in K, t in T, p in P))

        # non-existent vars
        @constraint(m, [c in C, j in DP, t in setdiff(Tendall, Tend[c])], dexp[c,j,t] == 0)

        @constraint(m, demandexpire[c in C, j in DP, t in Tend[c]], sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t) -
            sum(deltaQ[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t) <= dexp[c,j,t])
        @constraint(m, nonantecipativity[c in C, j in DP, t in T],
            sum(deltaQ[c,j,tau] for tau in 1:t) <= sum(demand[c,j,tau] for tau in 1:t))

        for c in C
            for j in DP
                for j2 in DP
                    for t in Tend[c]
                        if sum(demand[c,j,tau] for tau in t-threshold[c]:t) > ZERO && sum(demand[c,j2,tau] for tau in t-threshold[c]:t) > ZERO
                            @constraint(m, dexp[c,j,t]/sum(demand[c,j,tau] for tau in t-threshold[c]:t)
                                - dexp[c,j2,t]/sum(demand[c,j2,tau] for tau in t-threshold[c]:t) <= er[c,t])
                        end
                    end
                end
            end
        end

        # non-existent vars
        @constraint(m, [c in C, t in setdiff(Tendall, Tend[c])], er[c,t] == 0)
    end

    @constraint(m, supplycapacity[c in C, t in Tp], sum(deltaQ[c,j,tau] for j in DP, tau in t1:t)
        <= sum(supply[c,tau] for tau in 1:t)) # PS run: supply and demand are residual

    # non-existent vars
    @constraint(m, [c in C, i in N, k in K, t in Tp, p in P], Q[c,i,i,k,t,p] == 0)
    @constraint(m, [i in N, k in K, t in Tp, p in P], V[i,i,k,t,p] == 0)
    @constraint(m, [j in DP, k in K, t in Tp, p in P], V["theta",j,k,t,p] + V[j,"theta",k,t,p] == 0)
    @constraint(m, [i in N, j in N, k in K, t in Tp, i != j], V[i,j,k,t,1] == 0) # moment 1 is reserved to leave theta
    @constraint(m, sum(V["theta",o,k,t,p] for k in K, t in Tp, p in 2:nP) == 0)

    @constraint(m, [i in N, j in N, k in K, t in Tp, p in P; i != j], V[i,j,k,t,p] <= vehicle[k,t])
    @constraint(m, [k in K, t in Tp], V["theta",o,k,t,1] == sum(V[o,"theta",k,t,p] for p in 3:nP))
    @constraint(m, [k in K, t in Tp, p in 2:nP], sum(V[j,o,k,t,p-1] for j in ["theta";DP]) == sum(V[o,j,k,t,p] for j in ["theta";DP]))
    @constraint(m, [j in DP, k in K, t in Tp, p in P], sum(V[i,j,k,t,p] for i in filter(e -> e != j, N)) ==
        sum(V[j,i,k,t,p] for i in filter(e -> e != j, N)))
    @constraint(m, [k in K, t in Tp], sum(distance[i,j,t]/vspeed[vmodel[k]] * sum(V[i,j,k,t,p] for p in P)
        for i in N, j in N) + vloadtime[vmodel[k]] * sum(V[o,j,k,t,p] for j in DP, p in P) <= working_hours)
    @constraint(m, [i in N, j in N, k in K, t in Tp, p in P; i != j], sum(weight[c] * Q[c,i,j,k,t,p] for c in C) <= vload[vmodel[k]] * V[i,j,k,t,p])
    @constraint(m, sum(Q[c,j,o,k,t,p] for c in C, j in DP, k in K, t in Tp, p in P) == 0) # cut
    @constraint(m, [c in C, j in DP, t in Tp], sum(Q[c,i,j,k,t,p] - Q[c,j,i,k,t,p] for i in N, k in K, p in P) == deltaQ[c,j,t])
    @constraint(m, [c in C, j in DP, k in K, t in Tp, p in P], sum(Q[c,i,j,k,t,p] - Q[c,j,i,k,t,p] for i in N) >= 0) # sub-tour elimination constraints

    write_to_file(m, "./model" * instance * ".lp")

    return m
end

function proximitySearch(t1, t2, supplyps, demandhc, loadshc, renthc, fuelhc, vset, qdict)
    tstart = time()
    vsetps, qdictps = Set(), Dict()
    demandps = copy(demandhc)
    loadsps = copy(loadshc)
    fuelps = copy(fuelhc)
    rentps = Dict{Integer, Float64}()
    for t in T
        rentps[t] = sum(renthc[k,t] for k in K)
    end
    tlastimprovement = 0

    copySolution(C, No, K, t1-2:t1-1, vset, qdict, vsetps, qdictps)
    addLoad(t1-2, demandps, loadshc)
    addLoad(t1-1, demandps, loadshc)
    incumbent = sum(renthc[k,t] for k in K for t in t1-2:t1-1) + sum(fuelhc[t] for t in t1-2:t1-1) +
        sum(priority[c] * demandhc[c,j,t] for c in C for j in DP for t in 1:t1-1)

    for t in t1:t2
        copySolution(C, No, K, [t], vset, qdict, vsetps, qdictps)
        addLoad(t, demandps, loadshc)
        incumbent += sum(renthc[k,t] for k in K) + fuelhc[t] + sum(priority[c] * demandhc[c,j,t] for c in C for j in DP)

        # if there is unsatisfied demand, try 12% improvement
        # if sum(priority[c] * demandhc[c,j,tau] for c in C for j in DP for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:tend[c,t]) > minimum(values(vcost))
        #     incumbentaux = proximitySearch(C, DP, N, K, t-2, t, supplyps, demandps, loadsps, rentps, fuelps, vsetps, qdictps, incumbent, 0.88)
        #     if incumbentaux > 0
        #         incumbent = incumbentaux
        #         println("+++", incumbent)
        #         tlastimprovement = t
        #     end
        # end
        incumbentaux = proximitySearch(C, DP, N, K, t-2, t, supplyps, demandps, loadsps, rentps, fuelps, vsetps, qdictps, incumbent)
        if incumbentaux > 0
            incumbent = incumbentaux
            println("+++", incumbent)
            tlastimprovement = t
        elseif tlastimprovement == 0  # incumbent is an UB. taking period t-2 out to tighten the UB
            incumbent -= (sum(renthc[k,t-2] for k in K) + fuelhc[t-2] + sum(priority[c] * demandhc[c,j,t-2] for c in C for j in DP))
        end

        updateResidual(C, DP, t-2, supplyps, demandps, loadsps)
        # printRouting(C, No, K, T, P, vsetps, qdictps, vload, vmodel, demandps, 0)

        GC.gc()
        # flush(stdout)
        sleep(SLEEP)
    end

    updateResidual(C, DP, t2-1, supplyps, demandps, loadsps)
    updateResidual(C, DP, t2, supplyps, demandps, loadsps)

    pseq = getEquity(demandps)
    elapsedtimeps = time()-tstart

    inuse = 0
    for k in K
        if length(filter(e->e[3]==k && e[5]==1 && e[1]=="theta" && e[2]==o, vsetps)) > 0
            inuse += 1
        end
    end

    # if DEBUG
    #     printRouting(C, No, K, T, P, vsetps, qdictps, vload, vmodel, demandps, elapsedtimeps)
    # end

    return demandps, rentps, fuelps, sum(values(pseq)), Statistics.mean(values(pseq)), maximum(values(pseq)), elapsedtimeps, vsetps, qdictps, loadsps, inuse
end

function proximitySearch(C, DP, N, K, t1, t2, supply, demand, loads, rent, fuel, v, q, incumbent)
    Tp = T[t1:t2]
    Texp = []
    for c in C
        Texp = union(Texp,[tend[c,t1],tend[c,t2]])
    end

    m = createModel(C, DP, N, K, t1, t2, supply, demand, TIME_LIMIT)
    V = m[:V]
    Q = m[:Q]
    deltaQ = m[:deltaQ]
    dexp = m[:dexp]

    # We want to keep track of the original objective function:
    @expression(m, objExpr, sum(priority[c] * dexp[c,j,t] for c in C, j in DP, t in Texp)
        + sum(vcost[vmodel[k]] * V["theta",o,k,t,1] for k in K, t in Tp)
        + sum(vfuel[vmodel[k]] * distance[i,j,t] * V[i,j,k,t,p] for i in N, j in N, k in K, t in Tp, p in P))

    # Loop defining the iterations
    varOnes, varZeros = [], []
    iterations = 0
    while(true)
        iterations += 1
        # Separating variables:
        if iterations == 1
            varOnes, varZeros = separateZerosOnes(m, No, K, Tp, v, iterations == 1)
        end
        @show t1, iterations, size(varOnes), size(varZeros), incumbent

        # Redefining the objective function
        @objective(m, Min, sum(var for var in varZeros) + sum((1 - var) for var in varOnes))

        # Adding the Proximity Search constraint:
        @constraint(m, objExpr <= min(incumbent - 20, incumbent * IMPROVEMENT))
        set_name(@constraint(m, sum(var for var in varZeros) + sum((1 - var) for var in varOnes) <= MAXARCS), "maxarcs"*string(iterations))
        # if firstimprovement
        #     set_name(@constraint(m, sum(var for var in varZeros) + sum((1 - var) for var in varOnes) >= 4), "minarcs"*string(iterations))
        # end

        # if t2 == 4
        #     println(m)
            # printm = MathOptFormat.LP.Model()
            # MOI.copy_to(printm, JuMP.backend(m))
            # MOI.write_to_file(printm, "/Users/pcandida/Desktop/model.lp")
        # end

        # If we cannnot find a feasible solution, we stop!
        aux_time = @elapsed fo, rentt, fuelt, fairness, fairnessmax, gap, varOnes, varZeros = solveModel(m, C, DP, N, K, Tp, loads, v, q, 0)

        if fo < 0
            println(iterations, " ", aux_time)
            break
        end

        for t in Tp
            rent[t] = sum(vcost[vmodel[k]] * value(V["theta",o,k,t,1]) for k in K)
            fuel[t] = sum(vfuel[vmodel[k]] * distance[i,j,t] * ((i,j,k,t,p) in v ? 1 : 0)
                for i in N for j in filter(e -> e != i, N) for k in K for p in P)
        end

        # Getting value of the original objective function:
        mlp = subProblem(C, DP, N, K, t1, t2, supply, demand, getCapacity(v, N, K, Tp))
        aux = solveSubModel(mlp, C, DP, N, K, t1, t2, loads, dexp, Q, deltaQ, q)
        println(fo, " ", round(incumbent, digits=1), " ", round(value(objExpr), digits=1), " ", round(aux  + sum(rent[t] + fuel[t] for t in Tp), digits=1), " ", aux, " ", round(value(objExpr)-sum(rent[t] + fuel[t] for t in Tp), digits=1), " ", aux_time)
        incumbent = min(incumbent, value(objExpr), aux + sum(rent[t] + fuel[t] for t in Tp))

        delete(m, constraint_by_name(m, "maxarcs"*string(iterations)))
        # if firstimprovement
        #     delete(m, constraint_by_name(m, "minarcs"*string(iterations)))
        # end
    end

    return iterations > 1 ? incumbent : -1
end

function separateZerosOnes(m, N, K, T, v, first)
    varOnes, varZeros = [], []

    for i in N
        for j in N
            for k in K
                for p in P
                    for t in T
                        if (!first && value(m[:V][i,j,k,t,p]) > 0.5) || (first && (i,j,k,t,p) in v)
                            push!(varOnes, m[:V][i,j,k,t,p])
                            #println("V ", k, " ", t, " ", p, " (", i, " -> ", j, ")")
                        else
                            push!(varZeros, m[:V][i,j,k,t,p])
                        end
                    end
                end
            end
        end
    end

    return varOnes, varZeros
end

function getCapacity(v, N, K, T)
    vcapacity = Dict()
    for i in N
        for j in N
            for k in K
                for t in T
                    for p in P
                        vcapacity[i,j,k,t,p] = (i,j,k,t,p) in v ? vload[vmodel[k]] : 0
                    end
                end
            end
        end
    end
    return vcapacity
end

function solveModel(m, C, DP, N, K, T, loads, v, q, unsatisfiedDemand)
    varOnes = varZeros = rent = fuel = eq = eqmax = gap = fo = -1.0

    V = m[:V]
    Q = m[:Q]
    deltaQ = m[:deltaQ]

    optimize!(m)
    if has_values(m) && !isnan(objective_value(m))
        fo = objective_value(m)
        rent = sum(vcost[vmodel[k]] * value(V["theta",o,k,t,1]) for k in K for t in T)
        fuel = sum(vfuel[vmodel[k]] * distance[i,j,t] * value(V[i,j,k,t,p])
            for i in N for j in filter(e -> e != i, N) for k in K for t in T for p in P)
        gap = 0 # get_relative_mip_gap(m) in the next version of CPLEX.jl

        getRoutingvars(C, DP, N, K, T, V, Q, v, q)

        varOnes, varZeros = separateZerosOnes(m, No, K, T, v, false)

        if METHOD == EXACT
            eq = Statistics.mean(value(m[:er][c,t]) for c in C for t in Tend[c])
            eqmax = maximum(value(m[:er][c,t]) for c in C for t in Tend[c])
            for c in C
                for j in DP
                    for t in Tend[c]
                        unsatisfiedDemand[c,j,t] = value(m[:dexp][c,j,t])
                    end
                end
            end
        else
            for j in DP
                for c in C
                    for t in T
                        loads[c,j,t] = value(deltaQ[c,j,t])
                    end
                end
            end
            fuel += eliminateArcs(v, q, loads)
        end

    elseif METHOD == EXACT
        fo = objective_bound(m)
    end
    if DEBUG
        println(termination_status(m))
    end

    return fo, rent, fuel, eq, eqmax, gap, varOnes, varZeros
end

# PS and EXACT
function getRoutingvars(C, DP, N, K, T, V, Q, v, q)
    if length(v) > 0
        for i in ["theta";N]
            for j in ["theta";N]
                for k in K
                    for t in T
                        for p in P
                            if (i,j,k,t,p) in v
                                delete!(v, (i,j,k,t,p))
                                if i != "theta" && j != "theta" && j != o
                                    for c in C
                                        q[c,i,j,k,t,p] = 0
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    for i in ["theta";N]
        for j in ["theta";N]
            for k in K
                for t in T
                    for p in P
                        if value(V[i,j,k,t,p]) > 0.5
                            push!(v, (i,j,k,t,p))
                            # print("V ", k, " ", t, " ", p, " (", i, " -> ", j, ")")
                            if i != "theta" && j != "theta" && j != o
                                for c in C
                                    q[c,i,j,k,t,p] = value(Q[c,i,j,k,t,p])
                                    # @printf("\t%s = %.2f\t ", c, q[c,i,j,k,t,p])
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

function eliminateArcs(v, q, loads)
    fuel = 0
    for p in 2:nP
        for k in K
            for t in T
                for j in DP
                    if (o,j,k,t,p) in v && sum(loads[c,j,t] for c in C) < ZERO
                        for j2 in DP
                            if (j,j2,k,t,p) in v && pie[o,j2,t] != 0 && distance[o,j,t]+distance[j,j2,t] > distance[o,j2,t]
                                delete!(v, (o,j,k,t,p))
                                delete!(v, (j,j2,k,t,p))
                                delete!(v, (j2,j,k,t,p))
                                delete!(v, (j,o,k,t,p))
                                push!(v, (o,j2,k,t,p))
                                push!(v, (j2,o,k,t,p))
                                fuel += vfuel[vmodel[k]] * (distance[o,j2,t] + distance[j2,o,t] - distance[o,j,t] - distance[j,j2,t] - distance[j2,j,t] - distance[j,o,t])
                                for c in C
                                    q[c,o,j2,k,t,p] = q[c,o,j,k,t,p]
                                    q[c,o,j,k,t,p] = q[c,j,j2,k,t,p] = 0
                                end
                            end
                        end
                    end

                    for j2 in DP
                        if (j,j2,k,t,p) in v && (j2,j,k,t,p) in v && sum(q[c,j,j2,k,t,p] for c in C) < ZERO && sum(q[c,j2,j,k,t,p] for c in C) < ZERO
                            delete!(v, (j,j2,k,t,p))
                            delete!(v, (j2,j,k,t,p))
                            fuel -= vfuel[vmodel[k]] * (distance[j,j2,t] + distance[j2,j,t])
                        end
                    end
                end
            end
        end
    end
    return fuel
end

function calculateCosts(C, DP, N, K, T, P, v)
    rent = fuel = 0.0
    for k in K
        for t in T
            if ("theta",o,k,t,1) in v
                rent += vcost[vmodel[k]]
            end

            for p in P
                for i in N
                    for j in N
                        if (i,j,k,t,p) in v
                            fuel += vfuel[vmodel[k]] * distance[i,j,t]
                        end
                    end
                end
            end
        end
    end

    return rent, fuel
end

function subProblem(C, DP, N, K, t1, t2, supply, demand, vcapacity)
    Tp = T[t1:t2]
    Texp = []
    for c in C
        Texp = union(Texp,[tend[c,t1],tend[c,t2]])
    end

    m = Model(CPLEX.Optimizer)
    set_optimizer_attribute(m, "CPX_PARAM_MIPDISPLAY", 0)

    @variable(m, Q[C,N,N,K,Tp,P] >= 0)
    @variable(m, deltaQ[C,DP,Tp] >= 0)
    @variable(m, dexp[C,DP,Texp] >= 0)

    @objective(m, Min, sum(priority[c] * dexp[c,j,t] for c in C, j in DP, t in Texp))

    @constraint(m, [i in N, j in N, k in K, t in Tp, p in P; i != j], sum(weight[c] * Q[c,i,j,k,t,p] for c in C)
        <= vcapacity[i,j,k,t,p])

    @constraint(m, supplycapacity[c in C, t in Tp], sum(deltaQ[c,j,tau] for j in DP, tau in t1:t)
        <= sum(supply[c,tau] for tau in 1:t)) # PS run: supply and demand are residual

    @constraint(m, [c in C, i in N, k in K, t in Tp, p in P], Q[c,i,i,k,t,p] == 0) # non-existent var
    @constraint(m, sum(Q[c,j,o,k,t,p] for c in C, j in DP, k in K, t in Tp, p in P) == 0) # cut
    @constraint(m, [c in C, j in DP, k in K, t in Tp, p in P], sum(Q[c,i,j,k,t,p] - Q[c,j,i,k,t,p] for i in N) >= 0)
    @constraint(m, [c in C, j in DP, t in Tp], sum(Q[c,i,j,k,t,p] - Q[c,j,i,k,t,p] for i in N, k in K, p in P) == deltaQ[c,j,t])

    @constraint(m, demandexpire[c in C, j in DP, t in unique([tend[c,t1],tend[c,t2]])],
        sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:min(t,t2)) -
        sum(deltaQ[c,j,tau] for tau in max(div(t-1,threshold[c]+1)*(threshold[c]+1)+1,t1):min(t,t2)) <= dexp[c,j,t])
    @constraint(m, nonantecipativity[c in C, j in DP, t in setdiff(Tp, unique([tend[c,t1],tend[c,t2]]))],
        sum(deltaQ[c,j,tau] for tau in max(div(t-1,threshold[c]+1)*(threshold[c]+1)+1,t1):t)
        <= sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t))

    return m
end

function solveSubModel(m, C, DP, N, K, t1, t2, loads, dexp, Q, deltaQ, q)
    optimize!(m)
    if has_values(m) && !isnan(objective_value(m))
        dexplp = m[:dexp]
        Qlp = m[:Q]
        deltaQlp = m[:deltaQ]

        for c in C
            for j in DP
                for t in t1:t2
                    loads[c,j,t] = value(deltaQlp[c,j,t])
                    # set_start_value(deltaQ[c,j,t], loads[c,j,t])
                    # println(c, " ", j, " ", t, " ", loads[c,j,t])
                end

                # for t in unique([tend[c,t1],tend[c,t2]])
                #     set_start_value(dexp[c,j,t], value(dexplp[c,j,t]))
                # end
            end

            for j in N
                for i in N
                    for k in K
                        for t in t1:t2
                            for p in P
                                q[c,i,j,k,t,p] = value(Qlp[c,i,j,k,t,p])
                                # set_start_value(Q[c,i,j,k,t,p], q[c,i,j,k,t,p])
                                # println(q[c,i,j,k,t,p])
                            end
                        end
                    end
                end
            end
        end

        return objective_value(m)
    else
        return B
    end
end

# function fixSolutionRouting(C, N, K, T, P, V, Q, v, q)
#     if length(v) > 0
#         for k in K
#             for t in T
#                 for p in P
#                     for i in N
#                         for j in N
#                             aux = in((i,j,k,t,p), v)
#                             JuMP.fix(V[i,j,k,t,p], aux ? 1 : 0)
#                             # if aux print("V_", t, ", ", p, " (", i, " -> ", j, ")") end
#                             if i != "theta" && j != "theta"
#                                 for c in C
#                                     # if aux && in((c,i2,i,j,k,t,p), keys(q)) && q[c,i2,i,j,k,t,p] > 0
#                                     #     @show c, i, j, k, t, p, q[c,i2,i,j,k,t,p]
#                                     # end
#                                     JuMP.fix(Q[c,i,j,k,t,p], aux && in((c,i,j,k,t,p), keys(q)) ?
#                                         trunc(q[c,i,j,k,t,p],3) : 0)
#                                 end
#                             end
#                             # if aux println() end
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     flush(STDOUT)
# end
#
# function saveSolution(C, FN, N, K, T, P, v, q)
#     filename = string(filepath, instance, collaborative ? "c" : "n", "routing.txt")
#     open(filename, "w") do file
#         for i in N
#             for j in N
#                 for k in K
#                     for t in T
#                         for p in P
#                             aux = in((i,j,k,t,p), v)
#                             write(file, aux ? "1\n" : "0\n")
#                             if i != "theta" && j != "theta" && aux
#                                 for c in C
#                                     write(file, in((c,i,j,k,t,p), keys(q)) ? "$(q[c,i,j,k,t,p])\n" : "0\n")
#                                 end
#                             end
#                         end
#                     end
#                 end
#             end
#         end
#     end
# end
#
# function readSolution(C, FN, N, K, T, P)
#     filename = string(filepath, instance, collaborative ? "c" : "n", "routing.txt")
#     v = Set()
#     q = Dict()
#     open(filename, "r") do file
#         for i in N
#             for j in N
#                 for k in K
#                     for t in T
#                         for p in P
#                             aux = parse(UInt8, readline(file))
#                             if aux == 1
#                                 push!(v, (i,j,k,t,p))
#                                 if i != "theta" && j != "theta"
#                                     for c in C
#                                         q[c,i,j,k,t,p] = parse(Float64, readline(file))
#                                     end
#                                 end
#                             end
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     return v, q
# end
