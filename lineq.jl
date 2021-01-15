function linearEq(loads, vset, qdict, dexpcost, t1, t2)
    demandlin = copy(demand)
    supplylin = copy(supply)
    loadslin = copy(loads)
    vsetlin = copy(vset)
    qdictlin = copy(qdict)

    mlp = subProblemEq(supply, demand, getCapacity(vset), vsetlin, qdictlin, dexpcost, t1, t2)
    aux = solveSubModelEq(mlp, qdictlin, loadslin)

    for t in T
        updateResidual(C, DP, t, supplylin, demandlin, loadslin)
    end
    GC.gc()

    linfo = sum(priority[c] * demandlin[c,j,t] for c in C for j in DP for t in T)
    lineq = getEquity(demandlin)

    if DEBUG
        printRouting(C, No, K, T, P, vsetlin, qdictlin, vload, vmodel, demandlin, 0)
    end

    return linfo, sum(values(lineq)), qdictlin
end

function subProblemEq(supply, demand, vcapacity, v, q, dexpcost, t1, t2)
    m = Model(CPLEX.Optimizer)
    set_optimizer_attribute(m, "CPX_PARAM_MIPDISPLAY", 0)

    @variable(m, Q[C,N,N,K,T,P] >= 0)
    @variable(m, deltaQ[C,DP,T] >= 0)
    @variable(m, dexp[C,DP,Tendall] >= 0)
    @variable(m, 0 <= er[C,Tendall] <= 1)

    @objective(m, Min, sum(priority[c] * er[c,t] for c in C, t in Tend[c]))

    for c in C
        for t in Tend[c]
            t0 = div(t-1,threshold[c]+1)*(threshold[c]+1)+1
            for j in DP
                for j2 in DP
                    if sum(demand[c,j,tau] for tau in t0:t) > ZERO && sum(demand[c,j2,tau] for tau in t0:t) > ZERO
                        @constraint(m, dexp[c,j,t]/sum(demand[c,j,tau] for tau in t0:t)
                            - dexp[c,j2,t]/sum(demand[c,j2,tau] for tau in t0:t) <= er[c,t])
                    end
                end
            end
        end
    end

    # fix vars - previous optimisations
    for t in t1:t2
        for k in K
            for p in P
                for i in N
                    for j in DP
                        if in((i,j,k,t,p), v)
                            for c in C
                                @constraint(m, Q[c,i,j,k,t,p] == q[c,i,j,k,t,p])
                            end
                        end
                    end
                end
            end
        end
    end

    # non-existent vars
    @constraint(m, [c in C, t in setdiff(Tendall, Tend[c])], er[c,t] == 0)
    @constraint(m, [c in C, j in DP, t in setdiff(Tendall, Tend[c])], dexp[c,j,t] == 0)

    @constraint(m, sum(priority[c] * dexp[c,j,t] for c in C, j in DP, t in Tend[c]) <= dexpcost)
    # @constraint(m, sum(vcost[vmodel[k]] * V["theta",o,k,t,1] for k in K, t in T) +
    #     sum(vfuel[vmodel[k]] * pie[i,j,t] * distance[i,j] * V[i,j,k,t,p] for i in N for j in filter(e -> e != i, N) for k in K for t in T for p in P) <= transpcost)

    @constraint(m, [i in N, j in N, k in K, t in T, p in P; i != j], sum(weight[c] * Q[c,i,j,k,t,p] for c in C)
        <= vcapacity[i,j,k,t,p]+ZERO)

    @constraint(m, supplycapacity[c in C, t in T], sum(deltaQ[c,j,tau] for j in DP, tau in 1:t)
        <= sum(supply[c,tau] for tau in 1:t)+ZERO) # PS run: supply and demand are residual

    @constraint(m, [c in C, i in N, k in K, t in T, p in P], Q[c,i,i,k,t,p] == 0) # non-existent var
    @constraint(m, sum(Q[c,j,o,k,t,p] for c in C, j in DP, k in K, t in T, p in P) == 0) # cut
    @constraint(m, [c in C, j in DP, k in K, t in T, p in P], sum(Q[c,i,j,k,t,p] - Q[c,j,i,k,t,p] for i in N) >= 0)
    @constraint(m, [c in C, j in DP, t in T], sum(Q[c,i,j,k,t,p] - Q[c,j,i,k,t,p] for i in N, k in K, p in P) == deltaQ[c,j,t])

    @constraint(m, demandexpire[c in C, j in DP, t in Tend[c]], sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t) -
        sum(deltaQ[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t) == dexp[c,j,t])
    @constraint(m, nonantecipativity[c in C, j in DP, t in T],
        sum(deltaQ[c,j,tau] for tau in 1:t) <= sum(demand[c,j,tau] for tau in 1:t))

    return m
end

function solveSubModelEq(m, q, loads)
    optimize!(m)
    if has_values(m) && !isnan(objective_value(m))
        for i in N
            for j in N
                for k in K
                    for t in T
                        for p in P
                            for c in C
                                q[c,i,j,k,t,p] = value(m[:Q][c,i,j,k,t,p])
                            end
                        end
                    end
                end
            end
        end

        for j in DP
            for c in C
                for t in T
                    loads[c,j,t] = value(m[:deltaQ][c,j,t])
                end
            end
        end

        return sum(value(m[:er][c,t]) for c in C for t in Tend[c])#dexplp
    else
        return B
    end
end

function getCapacity(v)
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
