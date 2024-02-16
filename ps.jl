function setObjective(m)
    set_optimizer_attribute(m, "CPX_PARAM_MIPDISPLAY", 2)

    if EQCONSTR
        @objective(m, Min, sum(priority[c] * m[:dexp][c,j,t] for c in C, j in DP, t in Tend[c])
            + wfairness * sum(m[:er][c,t] for c in C, t in Tend[c])
            + wtransp * (sum(vcost[vmodel[k]] * m[:V]["theta",o,k,t,1] for k in K, t in T)
            + sum(vfuel[vmodel[k]] * distance[i,j,t] * m[:V][i,j,k,t,p] for i in N, j in N, k in K, t in T, p in P)))
    else
        @objective(m, Min, sum(priority[c] * m[:dexp][c,j,t] for c in C, j in DP, t in Tend[c])
            + sum(vcost[vmodel[k]] * m[:V]["theta",o,k,t,1] for k in K, t in T)
            + sum(vfuel[vmodel[k]] * distance[i,j,t] * m[:V][i,j,k,t,p] for i in N, j in N, k in K, t in T, p in P))
    end
end

function createModel(method)
    m = Model(CPLEX.Optimizer)

    @variable(m, V[["theta";N],["theta";N],K,T,P], Bin)
    @variable(m, Q[C,N,N,K,T,P] >= 0)
    @variable(m, deltaQ[C,DP,T] >= 0)
    @variable(m, dexp[C,DP,Tendall] >= 0)

    set_optimizer_attribute(m, "CPX_PARAM_THREADS", NTHREADS)
    set_optimizer_attribute(m, "CPX_PARAM_EPGAP", 0.001)
    set_optimizer_attribute(m, "CPX_PARAM_TILIM", method == PS ? MY_TIME_LIMIT : 3600)

    if method == PS
        set_optimizer_attribute(m, "CPX_PARAM_MIPDISPLAY", 0)
        set_optimizer_attribute(m, "CPXPARAM_MIP_Limits_Solutions", 1)
        set_optimizer_attribute(m, "CPXPARAM_MIP_Cuts_MIRCut", 1)
        set_optimizer_attribute(m, "CPXPARAM_MIP_Strategy_SubAlgorithm", 1)
        set_optimizer_attribute(m, "CPXPARAM_MIP_Strategy_HeuristicFreq", -1)
        set_optimizer_attribute(m, "CPXPARAM_MIP_Limits_CutPasses", -1)
    else#if method == FIXOPT
        set_optimizer_attribute(m, "CPX_PARAM_FLOWCOVERS", 1)
        set_optimizer_attribute(m, "CPX_PARAM_MIRCUTS", 1)
        set_optimizer_attribute(m, "CPX_PARAM_BTTOL", 0.1)
        set_optimizer_attribute(m, "CPX_PARAM_HEURFREQ", 100)
        set_optimizer_attribute(m, "CPX_PARAM_PROBE", -1)
    end

    # non-existent vars
    @constraint(m, [c in C, i in N, k in K, t in T, p in P], Q[c,i,i,k,t,p] == 0)
    @constraint(m, [i in No, k in K, t in T, p in P], V[i,i,k,t,p] == 0)
    @constraint(m, [j in DP, k in K, t in T, p in P], V["theta",j,k,t,p] + V[j,"theta",k,t,p] == 0)
    @constraint(m, [i in N, j in N, k in K, t in T, i != j], V[i,j,k,t,1] == 0) # moment 1 is reserved to leave theta
    @constraint(m, sum(V[o,"theta",k,t,1] for k in K, t in T) == 0)
    @constraint(m, sum(V["theta",o,k,t,p] for k in K, t in T, p in 2:nP) == 0)
    @constraint(m, [c in C, j in DP, t in setdiff(Tendall, Tend[c])], dexp[c,j,t] == 0)

    @constraint(m, [i in N, j in N, k in K, t in T, p in P; i != j], V[i,j,k,t,p] <= vehicle[k,t]*V["theta",o,k,t,1])
    @constraint(m, [k in K, t in T], V["theta",o,k,t,1] == sum(V[o,"theta",k,t,p] for p in 3:nP))
    @constraint(m, [k in K, t in T, p in 2:nP], sum(V[j,o,k,t,p-1] for j in ["theta";DP]) == sum(V[o,j,k,t,p] for j in ["theta";DP]))
    @constraint(m, [j in DP, k in K, t in T, p in P], sum(V[i,j,k,t,p] for i in filter(e -> e != j, N)) ==
        sum(V[j,i,k,t,p] for i in filter(e -> e != j, N)))
    @constraint(m, [k in K, t in T], sum(distance[i,j,t]/vspeed[vmodel[k]] * sum(V[i,j,k,t,p] for p in P)
        for i in N, j in N) + vloadtime[vmodel[k]] * sum(V[o,j,k,t,p] for j in DP, p in P) +
            sum(unloadtime * sum(V[i,j,k,t,p] for p in P) for i in N, j in DP) <= working_hours)

    @constraint(m, [i in N, j in N, k in K, t in T, p in P; i != j], sum(weight[c] * Q[c,i,j,k,t,p] for c in C) <= vload[vmodel[k]] * V[i,j,k,t,p])
    @constraint(m, sum(Q[c,j,o,k,t,p] for c in C, j in DP, k in K, t in T, p in P) == 0) # cut
    @constraint(m, [c in C, j in DP, t in T], sum(Q[c,i,j,k,t,p] - Q[c,j,i,k,t,p] for i in N, k in K, p in P) == deltaQ[c,j,t])
    @constraint(m, [c in C, j in DP, k in K, t in T, p in P], sum(Q[c,i,j,k,t,p] - Q[c,j,i,k,t,p] for i in N) >= 0) # sub-tour elimination constraints
    @constraint(m, supplycapacity[c in C, t in T], sum(deltaQ[c,j,tau] for j in DP, tau in 1:t)
        <= sum(supply[c,tau] for tau in 1:t))

    @constraint(m, demandexpire[c in C, j in DP, t in Tend[c]], sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t) -
        sum(deltaQ[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t) <= dexp[c,j,t])
    @constraint(m, nonantecipativity[c in C, j in DP, t in T],
        sum(deltaQ[c,j,tau] for tau in 1:t) <= sum(demand[c,j,tau] for tau in 1:t))

    if EQCONSTR
        @variable(m, 0 <= er[C,Tendall] <= 1)
        @constraint(m, [c in C, t in setdiff(Tendall, Tend[c])], er[c,t] == 0)
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
    end

    # creating LPs for exact model tuning (for PS tuning, objective function is different)
    # setObjective(m)
    # write_to_file(m, "./model/ps" * instance * ".lp")
    # exit(0)

    return m
end

function getRoutes(routes, routesch, routeindex, nt=1, routesvalue=Dict(), qr=Dict(), pastt=0)
    R = Dict()
    maxR = 0
    for t in T
        rt = length(routes[t])
        R[t] = collect(1:rt)
        maxR = max(maxR, rt)
    end
    allR = collect(1:maxR)

    m = Model(CPLEX.Optimizer)

    @variable(m, V2[K,T], Bin)
    @variable(m, V[K,allR,T], Bin)
    @variable(m, Q[C,K,allR,DP,T] >= 0)
    @variable(m, dexp[C,DP,Tendall] >= 0)

    set_optimizer_attribute(m, "CPX_PARAM_THREADS", NTHREADS)
    set_optimizer_attribute(m, "CPX_PARAM_EPGAP", 0.001)
    set_optimizer_attribute(m, "CPX_PARAM_TILIM", DYNAMIC_ANALYSIS ? 300 : nt*1800)
    # set_optimizer_attribute(m, "CPX_PARAM_MIPDISPLAY", 0)

    if EQCONSTR
        @variable(m, 0 <= er[C,Tendall] <= 1)
        @constraint(m, [c in C, t in setdiff(Tendall, Tend[c])], er[c,t] == 0)

        @objective(m, Min, sum(priority[c] * dexp[c,j,t] for c in C, j in DP, t in Tend[c])
            + wfairness * sum(er[c,t] for c in C, t in Tend[c])
            + wtransp * (sum(vcost[vmodel[k]] * V2[k,t] for k in K, t in T)
            + sum(vfuel[vmodel[k]] * getRouteDistance(routes[t][r],t) * V[k,r,t] for k in K, t in T, r in R[t])))
    else
        @objective(m, Min, sum(priority[c] * dexp[c,j,t] for c in C, j in DP, t in Tend[c])
            + sum(vcost[vmodel[k]] * V2[k,t] for k in K, t in T)
            + sum(vfuel[vmodel[k]] * getRouteDistance(routes[t][r],t) * V[k,r,t] for k in K, t in T, r in R[t]))
    end

    @constraint(m, [k in K, t in T, r in R[t]], V[k,r,t] <= V2[k,t])

    # non existing variables
    @constraint(m, [c in C, j in DP, t in setdiff(Tendall, Tend[c])], dexp[c,j,t] == 0)
    @constraint(m, [c in C, k in K, t in T, r in R[t], j in setdiff(DP, routes[t][r][2:end-1])], Q[c,k,r,j,t] == 0)
    for t in T
        diffR = setdiff(allR, R[t])
        if length(diffR) > 0
            @constraint(m, [k in K, r in diffR], V[k,r,t] == 0)
            @constraint(m, [c in C, k in K, j in DP, r in diffR], Q[c,k,r,j,t] == 0)
        end
    end

    @constraint(m, [k in K, t in T, r in R[t]], V[k,r,t] <= vehicle[k,t])
    @constraint(m, [k in K, t in T], sum((getRouteDistance(routes[t][r],t)/vspeed[vmodel[k]]
        + (length(routes[t][r])-2)*unloadtime + vloadtime[vmodel[k]]) * V[k,r,t] for r in R[t]) <= working_hours)

    @constraint(m, [k in K, t in T, r in R[t]], sum(weight[c] * Q[c,k,r,j,t] for c in C, j in routes[t][r][2:end-1]) <= vload[vmodel[k]] * V[k,r,t])
    @constraint(m, supplycapacity[c in C, t in T], sum(Q[c,k,r,j,tau] for k in K, j in DP, tau in 1:t, r in R[tau])
        <= sum(supply[c,tau] for tau in 1:t))

    @constraint(m, [c in C, j in DP, t in Tend[c]], sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t) -
        sum(Q[c,k,r,j,tau] for k in K, tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t, r in R[tau]) <= dexp[c,j,t])
    @constraint(m, nonantecipativity[c in C, j in DP, t in T],
        sum(Q[c,k,r,j,tau] for k in K, tau in 1:t, r in R[tau]) <= sum(demand[c,j,tau] for tau in 1:t))

    if EQCONSTR
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
    end

    # write_to_file(m, "./model/routes" * instance * ".lp")
    # exit(0)

    if routesch !== nothing
        for t in T
            for k in K
                for r in R[t]
                    set_start_value(V[k,r,t], 0)
                end
            end
        end
        for t in T
            for k in K
                for r in routeindex[k,t]:(routeindex[k,t]+length(routesch[k,t])-1)
                    # println(k, " ", t, " ", r, " ", routes[t][r])
                    set_start_value(V[k,r,t], 1)
                end
            end
        end
    end

    if pastt > 0
        for k in K
            for t in 1:pastt
                for r in R[t]
                    fix(V[k,r,t], routesvalue[k,r,t]; force = true)
                    if routesvalue[k,r,t] == 1
                        println(k, " ", r, " ", t)
                        for j in routes[t][r]
                            if j != o
                                print(j, " ")
                                for c in C
                                    if (c,k,r,j,t) in keys(qr)
                                        print(c, " ", qr[c,k,r,j,t], "\t")
                                        fix(Q[c,k,r,j,t], qr[c,k,r,j,t]; force = true)
                                    end
                                end
                                println()
                            end
                        end
                    end
                end
            end
        end
    end

    optimize!(m)
    vset, qdict, loads = Set(), Dict(), Dict()
    for c in C
        for j in DP
            for t in T
                loads[c,j,t] = 0
            end
        end
    end
    if has_values(m) && !isnan(objective_value(m))
    	# println("solveRoutes optimised")
    	# flush(stdout)
        eq = EQCONSTR ? [value(m[:er][c,t]) for c in C for t in Tend[c]] : getEquityValues(m[:dexp])
        unsatisfiedDemand = Dict()
        for c in C
            for j in DP
                for t in Tend[c]
                    unsatisfiedDemand[c,j,t] = value(m[:dexp][c,j,t])
                end
            end
        end

        for k in K
            for t in T
                p = 1
                for r in R[t]
                    if value(V[k,r,t]) > 0.9
                        if p == 1
                            push!(vset, ("theta",o,k,t,1))
                            p += 1
                        end

			for j in routes[t][r]
			    if j != o
                                for c in C
                                    qr[c,k,r,j,t] = 0
				end
                            end
                        end

                        u = o
                        j = ""
                        offsetindex = 1
                        while length(routes[t][r]) - offsetindex > 1
                            j = routes[t][r][end-offsetindex]
                            if sum(value(Q[c,k,r,j,t]) for c in C) > ZERO
                                push!(vset, (j,u,k,t,p))
                                for c in C
                                    qdict[c,j,u,k,t,p] = 0
                                end
                                break
                            else
                                offsetindex += 1
                            end
                        end
                        if j == ""
                            println("------------------ERROR no node served in route ", routes[t][r], " ------------------")
                        end
                        for i in reverse(routes[t][r][1:end-offsetindex-1])
                            if i == o || sum(value(Q[c,k,r,i,t]) for c in C) > ZERO
                                push!(vset, (i,j,k,t,p))
                                for c in C
                                    qr[c,k,r,j,t] = value(Q[c,k,r,j,t])
                                    qdict[c,i,j,k,t,p] = qr[c,k,r,j,t] + qdict[c,j,u,k,t,p]
                                    loads[c,j,t] += qr[c,k,r,j,t]
                                end
                                u = j
                                j = i
                            end
                        end
                        p += 1

                        routesvalue[k,r,t] = 1
                    else
                        routesvalue[k,r,t] = 0
                    end
                end
                if p > 1
                    push!(vset, (o,"theta",k,t,p))
                end
            end
        end

        return objective_value(m), sum(priorityo[c] * unsatisfiedDemand[c,j,t] for c in C for j in DP for t in Tend[c]),
            sum(vcost[vmodel[k]] * value(V2[k,t]) for k in K for t in T) +
            sum(vfuel[vmodel[k]] * getRouteDistance(routes[t][r],t) * value(V[k,r,t]) for k in K for t in T for r in R[t]) + remove_arcs(vset, qdict),
            sum(eq), vset, qdict, loads
    end
    # write_to_file(m, "route_problem"*instance*".lp")
    return -1, -1, -1, -1, vset, qdict, loads
end

# for each k unsed by constructive in t, set constraint forcing v_ijktp = 0 for all i,j,p
function removeUnusedVehicles(m, vset)
    V = m[:V]
    for k in K
        for t in T
            if !(("theta",o,k,t,1) in vset)
                @constraint(m, [i in N, j in N, p in P], V[i,j,k,t,p] == 0)
            end
        end
    end
end

function createModel()
    m = createModel(PS)
    # We want to keep track of the original objective function:
    if EQCONSTR
        @expression(m, objExpr, sum(priority[c] * m[:dexp][c,j,t] for c in C, j in DP, t in Tend[c])
            + wfairness * sum(m[:er][c,t] for c in C, t in Tend[c])
            + sum(vcost[vmodel[k]] * m[:V]["theta",o,k,t,1] for k in K, t in T)
            + sum(vfuel[vmodel[k]] * distance[i,j,t] * m[:V][i,j,k,t,p] for i in N, j in N, k in K, t in T, p in P))
    else
        @expression(m, objExpr, sum(priority[c] * m[:dexp][c,j,t] for c in C, j in DP, t in Tend[c])
            + sum(vcost[vmodel[k]] * m[:V]["theta",o,k,t,1] for k in K, t in T)
            + sum(vfuel[vmodel[k]] * distance[i,j,t] * m[:V][i,j,k,t,p] for i in N, j in N, k in K, t in T, p in P))
    end
    @constraint(m, PSconstr, objExpr <= B)
    return m
end

function proximitySearch(vset, incumbent, ntimes = B)
    tfirstimprovement, tlastimprovement = 1, nT
    auxtime = totaltime = 0
    better = false

    m = createModel()

    while tfirstimprovement > 0 && auxtime + totaltime < 3600 && ntimes > 0
        ntimes -= 1
        t1 = t2 = 0
        for t in tfirstimprovement:tlastimprovement
            auxtime = @elapsed incumbentaux = proximitySearch(m, t, vset, incumbent, totaltime)
            totaltime += auxtime
            if incumbentaux > 0
                better = true
                incumbent = incumbentaux
                if t1 == 0
                    t1 = t
                end
                t2 = t
            end
            # GC.gc()
            # sleep(SLEEP)
        end
        tfirstimprovement = t1 == 0 ? 0 : max(tend[C[1],t1] - WINDOW, 1)
        tlastimprovement  = t2 == 0 ? 0 : min(tend[C[1],t2] + WINDOW, nT)
    end

    if better
        return solveFixedModel(m, vset)
    end

    return -1*ones(8)
end

function solveFixedModel(m, vset, qdict=nothing, tfix=nothing)
    # delete(m, m[:PSconstr])
    set_normalized_rhs(m[:PSconstr], B)

    inuse = 0
    for k in K
        if length(filter(e->e[3]==k && e[5]==1 && e[1]=="theta" && e[2]==o, vset)) > 0
            inuse += 1
        end
    end

    fixVars(2*nT, m, vset, qdict, tfix)
    for i in No
        for j in No
            for p in P
                if isVar(i,j,p)
                    for k in K
                        for t in T
                            if !is_fixed(m[:V][i,j,k,t,p])
                                @show i,j,k,t,p
                            end
                        end
                    end
                end
            end
        end
    end
    println("solveFixedModel fixing done")
    flush(stdout)
    # printRouting(C, No, K, T, P, vset, qdict, vload, vmodel, demand, "PS fixed", 0, 0)

    setObjective(m)
    optimize!(m)
    if has_values(m) && !isnan(objective_value(m))
        eq = EQCONSTR ? [value(m[:er][c,t]) for c in C for t in Tend[c]] : getEquityValues(m[:dexp])
        qdict = getQ(m[:V], m[:Q])
        loads, unsatisfiedDemand = Dict(), Dict()
        for c in C
            for j in DP
                for t in Tend[c]
                    unsatisfiedDemand[c,j,t] = value(m[:dexp][c,j,t])
                end
                for t in T
                    loads[c,j,t] = value(m[:deltaQ][c,j,t])
                end
            end
        end

        return objective_value(m), sum(priority[c] * unsatisfiedDemand[c,j,t] for c in C for j in DP for t in Tend[c]),
            sum(vcost[vmodel[k]] * value(m[:V]["theta",o,k,t,1]) for k in K for t in T) +
            sum(vfuel[vmodel[k]] * distance[i,j,t] * value(m[:V][i,j,k,t,p]) for i in N for j in N for k in K for t in T for p in P) + remove_arcs(vset, qdict),
            sum(eq), inuse, unsatisfiedDemand, loads, qdict
    end
    write_to_file(m, "lp_problem"*instance*string(tfix===nothing ? 0 : tfix)*".lp")
    return -1*ones(8)
end

function proximitySearch(m, t, v, incumbent, totaltime, q=nothing)
    basic = q === nothing # non-basic == dynamic
    V = m[:V]
    Q = m[:Q]
    deltaQ = m[:deltaQ]
    dexp = m[:dexp]
    varOnes, varZeros = [], []
    iterations, aux_time = 0, 0.0
    rent, fuel, eqsum, unsatdemand, loads = -1, -1, nothing, nothing, nothing

    # Loop defining the iterations
    while aux_time + totaltime < 3600 - (nT-t+1)*MY_TIME_LIMIT/1.5
        iterations += 1
        fixVars(t, m, v, q, iterations)

        # Separating variables:
        varOnes, varZeros = separateZerosOnes(K, T, V, v)
        @show t, iterations, size(varOnes), size(varZeros), incumbent

        # Redefining the objective function
        @objective(m, Min, sum(var for var in varZeros) + sum((1 - var) for var in varOnes))

        # Adding the Proximity Search constraint:
        set_normalized_rhs(m[:PSconstr], min(incumbent - 20, incumbent * IMPROVEMENT))

        aux_time += @elapsed optimize!(m)
        # If we cannnot find a feasible solution, we stop!
        if !has_values(m) || isnan(objective_value(m))
            println(termination_status(m), " ", aux_time)
            delete(m, constraint_by_name(m, string(t)*"-"*string(iterations)))
            break
        end

        # println("PS: ", objective_value(m))
        updateArcs(V, v)
        aux, eqsum, unsatdemand, loads, q = sub_problem(getCapacity(v, N, K, T), t, q)

        rent = sum(vcost[vmodel[k]] * value(V["theta",o,k,t,1]) for t in T for k in K)
        fuel = sum(vfuel[vmodel[k]] * distance[i,j,t] * (value(V[i,j,k,t,p]) > 0.9 ? 1 : 0)
                for i in N for j in filter(e -> e != i, N) for k in K for t in T for p in P) +
                remove_arcs(v, q)

        # Getting value of the original objective function:
        vexpr = value(m[:objExpr])
        eqaux = EQCONSTR ? wfairness * eqsum : 0
        println(objective_value(m), " ", round(incumbent, digits=1), " ", round(vexpr, digits=1), " ",
            round(aux + rent + fuel + eqaux, digits=1), " ", aux_time)
        incumbent = min(incumbent, EQCONSTR ? incumbent : vexpr, aux + rent + fuel + eqaux)
        delete(m, constraint_by_name(m, string(t)*"-"*string(iterations)))
    end

    if basic
        return iterations > 1 ? incumbent : -1
    end
    return iterations > 1 ? (incumbent, eqsum, unsatdemand, loads, q, rent, fuel) : -1*ones(7)
end

function separateZerosOnes(K, T, V, v)
    varOnes, varZeros = [], []

    for i in No
        for j in No
            for p in P
                if isVar(i,j,p)
                    for k in K
                        for t in T
                            # if (!first && value(V[i,j,k,t,p]) > 0.9) || (first && (i,j,k,t,p) in v)
                            if (i,j,k,t,p) in v
                                push!(varOnes, V[i,j,k,t,p])
                                #println("V ", k, " ", t, " ", p, " (", i, " -> ", j, ")")
                            else
                                push!(varZeros, V[i,j,k,t,p])
                            end
                        end
                    end
                end
            end
        end
    end

    return varOnes, varZeros
end

function isVar(i,j,p)
    return !(i == j || (i == "theta" && j in DP) || (i in DP && j == "theta") || (i == "theta" && p > 1) || (p == 1 && (i != "theta" || j != o)))
end

function fixVars(topt, m, vset, qdict=nothing, tfix=topt)
    nodes = Set{String}()
    for i in No
        for j in No
            for p in P
                if isVar(i,j,p)
                    for k in K
                        for t in T
                            if t >= topt
                                if is_fixed(m[:V][i,j,k,t,p]) # freeing every arc at this period topt
                                    unfix(m[:V][i,j,k,t,p])
                                end
                                if t == topt && topt < nT && (i,j,k,t,p) in vset
                                    push!(nodes, i)
                                    push!(nodes, j)
                                end
                            else
                                fix(m[:V][i,j,k,t,p], (i,j,k,t,p) in vset ? 1 : 0; force = true)
                            end
                        end
                    end
                end
            end
        end
    end

    if topt >= 1 && topt <= nT
        nodes = setdiff(DP, nodes)
        set_name(@constraint(m, sum(m[:V][i,j,k,t,p] for i in nodes, j in nodes, k in K, t in topt:min(topt+WINDOW,nT), p in P) <= 4), string(topt)*"-"*string(tfix))
    end

    if qdict !== nothing
        fixQ(tfix, qdict, m[:Q])
    end
end

function fixQ(tnow, qdict, Q)
    for i in N
        for j in N
            for p in P
                if isVar(i,j,p)
                    for k in K
                        for t in 1:tnow-1
                            for c in C
                                fix(Q[c,i,j,k,t,p], (c,i,j,k,t,p) in keys(qdict) ? qdict[c,i,j,k,t,p] : 0; force = true)
                                # if (c,i,j,k,t,p) in keys(qdict) && qdict[c,i,j,k,t,p] > ZERO @show (c,i,j,k,t,p), qdict[c,i,j,k,t,p] end
                            end
                        end
                    end
                end
            end
        end
    end
end

function updateArcs(V, v)
    empty!(v)

    for k in K
        for t in T
            for p in P
                for i in No
                    for j in filter(e -> e != i, No)
                        if value(V[i,j,k,t,p]) > 0.9
                            push!(v, (i,j,k,t,p))
                        end
                    end
                end
            end
        end
    end
end

function getQ(V, Q)
    q = Dict()

    for i in N
        for j in filter(e -> e != i, N)
            for k in K
                for t in T
                    for p in P
                        if value(V[i,j,k,t,p]) > 0.9
                            for c in C
                                q[c,i,j,k,t,p] = value(Q[c,i,j,k,t,p])
                            end
                        end
                    end
                end
            end
        end
    end

    return q
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

function sub_problem(vcapacity, tnow, qdict)
    m = Model(CPLEX.Optimizer)
    set_optimizer_attribute(m, "CPX_PARAM_MIPDISPLAY", 0)

    @variable(m, Q[C,N,N,K,T,P] >= 0)
    @variable(m, deltaQ[C,DP,T] >= 0)
    @variable(m, dexp[C,DP,Tendall] >= 0)

    @objective(m, Min, sum(priority[c] * dexp[c,j,t] for c in C, j in DP, t in Tend[c]))

    @constraint(m, [i in N, j in N, k in K, t in T, p in P; i != j], sum(weight[c] * Q[c,i,j,k,t,p] for c in C)
        <= vcapacity[i,j,k,t,p])

    @constraint(m, supplycapacity[c in C, t in T], sum(deltaQ[c,j,tau] for j in DP, tau in 1:t)
        <= sum(supply[c,tau] for tau in 1:t))

    @constraint(m, [c in C, i in N, k in K, t in T, p in P], Q[c,i,i,k,t,p] == 0) # non-existent var
    @constraint(m, sum(Q[c,j,o,k,t,p] for c in C, j in DP, k in K, t in T, p in P) == 0) # cut
    @constraint(m, [c in C, j in DP, k in K, t in T, p in P], sum(Q[c,i,j,k,t,p] - Q[c,j,i,k,t,p] for i in N) >= 0)
    @constraint(m, [c in C, j in DP, t in T], sum(Q[c,i,j,k,t,p] - Q[c,j,i,k,t,p] for i in N, k in K, p in P) == deltaQ[c,j,t])

    @constraint(m, demandexpire[c in C, j in DP, t in Tend[c]], sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t) -
        sum(deltaQ[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t) <= dexp[c,j,t])
    @constraint(m, nonantecipativity[c in C, j in DP, t in T],
        sum(deltaQ[c,j,tau] for tau in 1:t) <= sum(demand[c,j,tau] for tau in 1:t))

    if qdict !== nothing
        fixQ(tnow, qdict, Q)
    else
        qdict = Dict()
    end

    optimize!(m)
    if has_values(m) && !isnan(objective_value(m))
        eq = getEquityValues(dexp)

        loads, unsatisfiedDemand = Dict(), Dict()
        for c in C
            for j in DP
                for t in Tend[c]
                    unsatisfiedDemand[c,j,t] = value(m[:dexp][c,j,t])
                end
                for t in T
                    loads[c,j,t] = value(m[:deltaQ][c,j,t])
                end
            end
        end

        for c in C
            for i in N
                for j in N
                    for k in K
                        for t in T
                            for p in P
                                qdict[c,i,j,k,t,p] = value(Q[c,i,j,k,t,p])
                            end
                        end
                    end
                end
            end
        end

        return objective_value(m), sum(eq), unsatisfiedDemand, loads, qdict #Statistics.mean(eq), maximum(eq),
    end
    write_to_file(m, "sub_problem"*instance*string(tnow)*".lp")
    return B, -1, nothing, nothing, nothing
end

function remove_arcs(v, q)
    fuel = 0
    for k in K
        for t in T
            for p in 2:nP
                i = o
                route = []
                arcs = Set()
                push!(route, o)

                added = true
                while added
                    added = false
                    for j in DP
                        if (i,j,k,t,p) in v && !((i,j) in arcs)
                            push!(route, j)
                            push!(arcs, (i,j))
                            i = j
                            added = true
                            break
                        end
                    end
                end

                if length(route) > 2
                    push!(route, o)
                    # println(route)
                    # printRouting(C, No, [k], [t], P, v, q, vload, vmodel, demand, "PS", 0, 0)
                    i = popfirst!(route)
                    j = popfirst!(route)
                    added = true

                    while length(route) > 0 && j != o
                        j2 = popfirst!(route)
                        same = true
                        for c in C
                            if q[c,i,j,k,t,p] != q[c,j,j2,k,t,p]
                                same = false
                                break
                            end
                        end
                        if same && distance[i,j2,t] < distance[i,j,t] + distance[j,j2,t]
                            delete!(v, (i,j,k,t,p))
                            delete!(v, (j,j2,k,t,p))
                            if i != j2
                                push!(v, (i,j2,k,t,p))
                            end
                            fuel += vfuel[vmodel[k]] * (distance[i,j2,t] - distance[i,j,t] - distance[j,j2,t])
                            # print("\n\n\n\n", k, " ", t, " ", p, " ", i, " ", j, " ", j2, " ", distance[i,j2,t], " ", distance[i,j,t], " ", distance[j,j2,t], " ")
                            for c in C
                                # print(q[c,i,j,k,t,p], " ", q[c,j,j2,k,t,p], " ")
                                q[c,i,j2,k,t,p] = q[c,j,j2,k,t,p]
                                q[c,i,j,k,t,p] = q[c,j,j2,k,t,p] = 0
                            end
                            # println("\n\n\n")
                        else
                            i = j
                        end
                        j = j2
                    end
                    # printRouting(C, No, [k], T, P, v, q, vload, vmodel, demand, "PS", 0, 0)
                end
                if !added
                    break
                end
            end
        end
    end
    # println("saindo\n\n\n")
    return fuel
end

function getEquityValues(dexp)
    eq, aux = [], []
    for c in C
        for t in Tend[c]
            for j in DP
                if sum(demand[c,j,tau] for tau in t-threshold[c]:t) > 1
                    for j2 in DP
                        if j2 != j && sum(demand[c,j2,tau] for tau in t-threshold[c]:t) > 1
                            push!(aux, value(dexp[c,j,t])/sum(demand[c,j,tau] for tau in t-threshold[c]:t)
                                - value(dexp[c,j2,t])/sum(demand[c,j2,tau] for tau in t-threshold[c]:t))
                        end
                    end
                end
            end
            push!(eq, maximum(aux))
            empty!(aux)
        end
    end
    return eq
    # return [maximum(value(dexp[c,j,t])/sum(demand[c,j,tau] for tau in t-threshold[c]:t)
    #     - value(dexp[c,j2,t])/sum(demand[c,j2,tau] for tau in t-threshold[c]:t) for j in DP for j2 in DP)
    #     for c in C for t in Tend[c]] # NaN when one j has no demand
end

function fixopt(vset, qdict)
    m = createModel(FIXOPT)
    setObjective(m)

    tfirstimprovement = WINDOW
    tlastimprovement = nT
    totaltime = 0
    while tfirstimprovement > 0 && totaltime < 3000
        t1 = t2 = 0
        for t in tfirstimprovement:tlastimprovement
            fixVars(t-WINDOW+1:t, m, vset)
            # if t % 3 == 0
            #     write_to_file(m, "model"*instance*string(t)*".lp")
            # end
            totaltime += @elapsed optimize!(m)
            if has_values(m) && !isnan(objective_value(m))
                getRoutingvars(C, DP, N, K, t-WINDOW+1:t, m[:V], m[:Q], vset, qdict) # update vset
                if t1 == 0
                    t1 = t-WINDOW+1
                end
                t2 = t
            end
        end
        tfirstimprovement = t1 == 0 ? 0 : max(t1, WINDOW)
        tlastimprovement = min(t2 + WINDOW - 1, nT)
    end
    fixVars(T, vset, m[:V])
    optimize!(m)
    return objective_value(m), objective_bound(m), Statistics.mean(value(m[:er][c,t]) for c in C for t in Tend[c])
end

function updateV(relax, t, m)
    for i in N#o
        for j in N#o
            if i != j
                for k in K
                    # if i != "theta" && j != "theta"
                        for p in P
                            relax ? unset_binary(m[:V][i,j,k,t,p]) : @constraint(m, m[:V][i,j,k,t,p] in MOI.ZeroOne())
                        end
                    # end
                end
            end
        end
    end
end

function relax_fix()
    vset = Set()
    m = createModel(RELAXFIX)
    setObjective(m)

    for t in 2:nT
        updateV(true, t, m)
    end

    for t in T
        optimize!(m)
        println(objective_value(m))

        if t != nT
            updateArcs(m[:V], vset)
            fixVars(t+1:nT, m, vset)

            updateV(false, t+1, m)
        end
    end

    println(objective_value(m))
end
