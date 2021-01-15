include("functions.jl")
include("input.jl")
include("data.jl")
include("constructive.jl")
include("ps.jl")
include("lineq.jl")
include("debug.jl")

hcrent, hcfuel, hcfo, hceqsum, hceqmean, hceqmax, hcelapsedtime = 0, 0, 0, 0, 0, 0, 0 #Array{Float64}(undef, 1), Array{Float64}(undef, 1), Array{Float64}(undef, 1), Array{Float64}(undef, 1), Array{Float64}(undef, 1), Array{Float64}(undef, 1), Array{Float64}(undef, 1)
psrent, psfuel, psfo, pseqsum, pseqmean, pseqmax, pselapsedtime = 0, 0, 0, 0, 0, 0, 0 #Array{Float64}(undef, 1), Array{Float64}(undef, 1), Array{Float64}(undef, 1), Array{Float64}(undef, 1), Array{Float64}(undef, 1), Array{Float64}(undef, 1), Array{Float64}(undef, 1)
hcset, psset, hcdict, psdict, hcloads, psloads = Set(), Set(), Dict(), Dict(), Dict(), Dict()
updated = Dict()


for t in 1:nT-2
    fixinterval = t != nT-2 ? [t] : [t:nT;]
    global demandhc
    global demandps

    if t == 1
        supplyhc, demandhc = copy(supply), copy(demand)
    else
        getDataUpdate(t, updated)

        for t in T
            Printf.@printf("-----%d------\n", t)
            for j in DP
                for c in C
                    if demand[c,j,t] > ZERO
                        Printf.@printf("%s %s %.1f\t", c, j, demand[c,j,t])
                    end
                end
                println()
            end
        end

        supplyhc, demandhc = getResidual(t, psloads)

        for t in T
            Printf.@printf("-----%d------\n", t)
            for j in DP
                for c in C
                    if demand[c,j,t] > ZERO
                        Printf.@printf("%s %s %.1f\t", c, j, demandhc[c,j,t])
                    end
                end
                println()
            end
        end
    end
    supplyps = copy(supplyhc)

    demandhc, loads, rent, fuel, hceqsum, hceqmean, hceqmax, vset, qdict, hcelapsedtime = createSolution("onlybig", t, nT, demandhc, supplyhc)
    r, f = updateValues(t, t == nT-2 ? nT : t, rent, fuel, loads, hcloads)
    global hcrent += r
    global hcfuel += f
    copySolution(C, No, K, fixinterval, vset, qdict, hcset, hcdict)
    # if DEBUG
    #     printRouting(C, No, K, fixinterval, P, hcset, hcdict, vload, vmodel, demandhc, hcelapsedtime)
    # end

    demandps, rent, fuel, pseqsum, pseqmean, pseqmax, aux, vset, qdict, loads, carsused = proximitySearch(t+2, t+2, supplyps, demandhc, loads, rent, fuel, hcset, hcdict)
    global pselapsedtime += aux
    r, f = updateValues(t, t, rent, fuel, loads, psloads)
    global psrent += r
    global psfuel += f
    copySolution(C, No, K, fixinterval, vset, qdict, psset, psdict)
end

# --------------------------------------------------------------------------- #
hcfo = sum(priority[c] * demandhc[c,j,t] for c in C for j in DP for t in T)
psfo = sum(priority[c] * demandps[c,j,t] for c in C for j in DP for t in T)

if DEBUG
    printRouting(C, No, K, T, P, psset, psdict, vload, vmodel, demandps, pselapsedtime)
end

print("\n", hcrent, " ", hcfuel, " ", hcfo, " ", hceqsum, " ", hceqmax, " ", hcelapsedtime)
print(" ", psrent, " ", psfuel, " ", psfo, " ", pseqsum, " ", pseqmax, " ", pselapsedtime)
# print(" ", linfo, " ", lineqsum)
