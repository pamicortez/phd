function printRouting(C, No, K, T, P, V, Q, vload, vmodel, unsatisfiedDemand, method, fo, fot=nothing, topt=nothing)
    for k in K
        println("..................................Loads for vehicle ", k, " (", vload[vmodel[k]]/1000, " ton)..................................")
        for t in T
            println("..............", t, "..............")
            for p in P
                for i in No
                    for j in No
                        if (i,j,k,t,p) in V
                            print("V_", p, " (", i, " -> ", j, ")")
                            if i != "theta" && j != "theta" && j != o
                                Printf.@printf("%d, %.1f", distance[i,j,t], distance[i,j,t]/vspeed[vmodel[k]])
                                for c in C
                                    if ((c,i,j,k,t,p) in keys(Q))
                                        # Printf.@printf("\t erro %s, %s, %s, %d, %d", i,j,k,t,p)
                                    # else
                                        Printf.@printf("\t%.1f", Q[c,i,j,k,t,p])
                                    end
                                end
                            end
                            println()
                        end
                    end
                end
            end
        end
    end

    println((topt !== nothing ? topt : ""), ".....................", (fot !== nothing ? fot : ""), " ", method, " ", fo, ".........................")
    flush(stdout)
end

function print_loads(loads)
    # n = maximum([length(j) for j in DP])
    Printf.@printf("%60s %30s %18s", C[1], C[2], C[3])
    println()
    for j in DP
        Printf.@printf("%36s", j)
        for c in C
            ini = 1
            for t in Tend[c]
                Printf.@printf("%12.2f", sum(loads[c,j,t] for tau in ini:t))
                ini = t+1
            end
        end
        println()
    end
end
