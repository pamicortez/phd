function printRouting(C, No, K, T, P, V, Q, vload, vmodel, unsatisfiedDemand, elapsedTime)
    for k in K
        println("..................................Loads for vehicle ", k, " (", vload[vmodel[k]]/1000, " ton)..................................")
        for t in T
            println("..............", t, "..............")
            for p in 1:nP+10
                for i in No
                    for j in No
                        if (i,j,k,t,p) in V
                            print("V_", p, " (", i, " -> ", j, ")")
                            if i != "theta" && j != "theta" && j != o
                                Printf.@printf("%d, %.1f", distance[i,j,t], distance[i,j,t]/vspeed[vmodel[k]])
                                for c in C
                                    Printf.@printf("\t\"%s\" = %.1f\t ", c, Q[c,i,j,k,t,p])
                                end
                            end
                            println()
                        end
                    end
                end
            end
        end
    end

    # D
    println()
    for j in DP
        Printf.@printf("\"%s\"\t", j)
    end
    println()
    for c in C
        for j in DP
            Printf.@printf("%.1f\t", sum(unsatisfiedDemand[c,j,t] for t in T))
        end
        println()
    end

    println(".....................Routing Problem running time ", elapsedTime, " s.........................")
    flush(stdout)
end
