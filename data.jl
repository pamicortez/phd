# update residual supply and demand
function updateResidual(C, DP, t, residualsupply, residualdemand, loads)
    updateResidualSupply(C, DP, t, residualsupply, loads)
    updateResidualDemand(C, DP, t, 1, residualdemand, loads)
end

function updateResidualSupply(C, DP, t, residualsupply, loads)
    for c in C
        for j in DP
            auxs = loads[c,j,t]
            if auxs > 0
                for t2 in t:-1:1
                    if auxs < ZERO
                        break
                    end
                    if residualsupply[c,t2] >= auxs
                        residualsupply[c,t2] -= auxs #trunc(Int, auxs + ZERO)
                        if residualsupply[c,t2] > -ZERO && residualsupply[c,t2] < ZERO
                            residualsupply[c,t2] = 0
                        end
                        auxs = 0
                    else#if residualsupply[c,t2] > ZERO
                        auxs -= residualsupply[c,t2]
                        residualsupply[c,t2] = 0
                    end
                end
            end
        end
    end
end

function updateResidualDemand(C, DP, t, direction, residualdemand, loads)
    for c in C
        for j in DP
            auxd = loads[c,j,t]
            if auxd > 0
                ti = direction == 1 ? div(t-1,threshold[c]+1)*(threshold[c]+1)+1 : t
                tf = direction == 1 ? t : div(t-1,threshold[c]+1)*(threshold[c]+1)+1
                for t2 in ti:direction:tf
                    if auxd < ZERO
                        break
                    end
                    if residualdemand[c,j,t2] >= auxd
                        residualdemand[c,j,t2] -= auxd #trunc(Int, auxd + ZERO)
                        if residualdemand[c,j,t2] > -ZERO && residualdemand[c,j,t2] < ZERO
                            residualdemand[c,j,t2] = 0
                        end
                        auxd = 0
                    else#if residualdemand[c,j,t2] > ZERO
                        auxd -= residualdemand[c,j,t2]
                        residualdemand[c,j,t2] = 0
                    end
                    # println(c, " ", j, " ", t2, " ", residualdemand[c,j,t2])
                end
            end
        end
    end
end

function getResidual(t, loads)
    residualsupply = copy(supply)
    residualdemand = copy(demand)

    for taux in 1:t-1
        updateResidual(C, DP, taux, residualsupply, residualdemand, loads)
    end

    return residualsupply, residualdemand
end

# threshold[c] >= subperiod
function initialDemand(subperiod, demandps, loadshc)
    for c in C
        for j in DP
            t = subperiod
            while t < nT
                suml = sum(loadshc[c,j,t2] for t2 in t:tend[c,t])
                t2 = tend[c,t]
                while suml > ZERO
                    if t2 < 1
                        println("ERRO!!!")
                    end
                    if demandps[c,j,t2] >= suml
                        demandps[c,j,t2] -= suml
                        suml = 0
                    else demandps[c,j,t2] > ZERO
                        suml -= demandps[c,j,t2]
                        demandps[c,j,t2] = 0
                    end
                    t2 -= 1
                end
                # println(c, " ", j, " ", sum(demandps[c,j,t] for t in 1:tend[c,1]))
                t = tend[c,t] + 1
            end
        end
    end
end

function addLoad(t, demandps, loadshc)
    for c in C
        for j in DP
            suml = loadshc[c,j,t]
            # print(c, " ", j, " ", suml)
            t2 = t
            while suml > ZERO
                if t2 < div(t-1,threshold[c]+1)*(threshold[c]+1)+1
                    println("ERRO!!!")
                end
                if demand[c,j,t2] >= suml
                    demandps[c,j,t2] += suml
                    suml = 0
                else demand[c,j,t2] > ZERO
                    suml -= demand[c,j,t2]
                    demandps[c,j,t2] += demand[c,j,t2]
                end
                t2 -= 1
            end
            # println(" ", sum(demandps[c,j,t] for t in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:tend[c,t]))
        end
    end
end

function preprocessing(t, residualsupply, residualdemand)
    sday = Dict{Any,Float64}()
    dday = Dict{Any,Float64}()

    # supply in t
    for c in C
        sday[c] = 0
    end
    for c in C
        sday[c] = sum(residualsupply[c,tau] for tau in 1:t)
        if sday[c] > -ZERO && sday[c] < ZERO
            sday[c] = 0
        end
        # println(c, " ", sday[c])
        if sday[c] < 0
            println("negative supply!!!!!!!!!!! ", c, " ", sday[c])
        end
    end

    # demand in t
    for c in C
        for j in DP
            dday[c,j] = sum(residualdemand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:t)
            if dday[c,j] > -ZERO && dday[c,j] < ZERO
                dday[c,j] = 0
            end
            if dday[c,j] < 0
                println("negative demand!!!!!!!!!!! ", c, " ", j, " ", dday[c,j])
            end
        end
    end

    return sday, dday
end

# CONSTRUCTIVE and CG
function getRoutingvars(m, k, t, v, q)
    for p in P
        for i in No
            for j in No
                if getvalue(m[:V][i,j,p]) > ZERO
                    push!(v, (i,j,k,t,p+10))
                    # if DEBUG print("V ", k, " ", t, " ", p, " (", i, " -> ", j, ")") end
                    if i != "theta" && j != "theta"
                        for c in C
                            q[c,i,j,k,t,p+10] = getvalue(m[:Q][c,i,j,p])
                            # if DEBUG @printf("\t%s = %.2f\t ", c, q[c,i,j,k,t,p]) end
                        end
                    end
                    # if DEBUG println() end
                end
            end
        end
    end
    return v, q
end

function getPrio(t, sday, dday, loads)
    prio = Dict{Any,Float64}()
    for j in DP
        for c in C
            if sday[c] > ZERO && dday[c,j] > ZERO
                prio[c,j] = priority[c] * (1 - sum(loads[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:tend[c,t])/
                    sum(demand[c,j,tau] for tau in div(t-1,threshold[c]+1)*(threshold[c]+1)+1:tend[c,t]))

                # if t == tend[c,t] # last day to deliver this order
                #     prio[c,j] *= 2
                # end
            else
                prio[c,j] = 0
            end
            # println(c, " ", j, " ", round(prio[c,j], digits=5))
        end
    end
    # flush(stdout)
    return prio
end

function getDPsortedByPrio(prio, dday)
    aux = copy(DP)
    i = 1
    prioaux = Array{Float64}(undef, nDP)
    for j in aux
        prioaux[i] = sum(dday[c,j] for c in C)
        i += 1
    end
    aux = aux[sort([1:nDP;], by=i->(prioaux[i]), rev=true)]

    i = 1
    for j in aux
        prioaux[i] = sum(prio[c,j] for c in C)
        i += 1
    end
    return aux[sort([1:nDP;], by=i->(prioaux[i]), rev=true)]
end

function getEquity(demand1)
    eq = Dict()
    for c in C
        for t in Tend[c]
            t1 = div(t-1,threshold[c]+1)*(threshold[c]+1)+1
            eq[c,t] = 0
            for j in DP
                for j2 in DP
                    if sum(demand[c,j,tau] for tau in t1:t) > ZERO && sum(demand[c,j2,tau] for tau in t1:t) > ZERO
                        eq[c,t] = max(eq[c,t], sum(demand1[c,j,tau] for tau in t1:t)/sum(demand[c,j,tau] for tau in t1:t)
                            - sum(demand1[c,j2,tau] for tau in t1:t)/sum(demand[c,j2,tau] for tau in t1:t))
                        # println(sum(demand1[c,j,tau] for tau in t1:t), " ", sum(demand[c,j,tau] for tau in t1:t), " ",
                        #    sum(demand1[c,j2,tau] for tau in t1:t), " ", sum(demand[c,j2,tau] for tau in t1:t))
                    end
                end
            end
        end
    end

    return eq #Statistics.mean(values(eq)), maximum(values(eq))
end

function copyAndFix(t1, t2, vset2, vset1, qdict2, qdict1, rent2, rent1, fuel2, fuel1, loads2, loads1, availableTime)
    for c in C
        for j in DP
            for t in t1:t2
                loads2[c,j,t] = loads1[c,j,t]
                for k in K
                    availableTime[k,t] = 0
                    rent2[k,t] = rent1[k,t]
                end
                fuel2[t] = fuel1[t]
                # tem que copiar todo vset tb, pois é usado em PS::separateZerosOnes
            end
        end
    end
    for k in K
        for t in t1:t2
            for p in P
                for i in No
                    for j in No
                        if (i,j,k,t,p) in vset1
                            push!(vset2, (i,j,k,t,p))
                            if i != "theta" && j != "theta" && j != o
                                for c in C
                                    qdict2[c,i,j,k,t,p] = qdict1[c,i,j,k,t,p]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return vset2, qdict2, rent2, fuel2, loads2
end

function setVars(loads, fuel, rent, availableTime, lastp)
    for c in C
        for j in DP
            for t in T
                loads[c,j,t] = 0
            end
        end
    end
    for t in T
        fuel[t] = 0
        for k in K
            rent[k,t] = 0
            availableTime[k,t] = working_hours
            lastp[k,t] = 0
        end
    end
end

# change reverberates to next periods
function updateDemand(C, DP, t, sign, demand, loads)
    for c in C
        for j in DP
            for t2 in tend[c,t]+1:threshold[c]+1:nT
                demand[c,j,t2] += sign*loads[c,j,t]
            end
        end
    end
end

function updateDemand(t, j, r, howmuch)
    virtualload = Dict()
    for c in C
        virtualload[c,j,t] = c == "cleaning kit" ? round(Int, howmuch/4) : howmuch
    end
    if r < 0
        updateResidualDemand(C, [j], t, -1, demand, virtualload)
    else
        for c in C
            demand[c,j,t] += virtualload[c,j,t]
        end
    end
    updateDemand(C, [j], t, r, demand, virtualload)
    # println(t, " ", j, " ", r, " ", howmuch)
end

function updateSupply(t, totalincrease)
    for c in C
        supply[c,t] += c == "cleaning kit" ? round(Int, totalincrease/4) : totalincrease
        for t2 in tend[c,t]+1:threshold[c]+1:nT
            supply[c,t2] += c == "cleaning kit" ? round(Int, totalincrease/4) : totalincrease
        end
    end
    # println(t, " updateSupply ", totalincrease)
end

function getDataUpdate(t, updated)
    filename = string(filepath, instance, "day", t, ".txt")
    if isfile(filename)
        local aux = readdlm(filename, '$')
        for i in 1:trunc(Int, length(aux)/3)
            if aux[i,2] == ""
                break
            end
            j = aux[i,1]
            r = aux[i,2]
            howmuch = aux[i,3]
            # println("demand ", demand["cleaning kit",j,1])
            updateDemand(t, j, r, howmuch)
            # println("demand2 ", demand["cleaning kit",j,1])
        end

        # println("supply ", supply["cleaning kit",1])
        updateSupply(t, length(aux) > 1 ? aux[trunc(Int, length(aux)/3),1] : aux[1])
        # println("supply2 ", supply["cleaning kit",1])
    else
        if isempty(updated) || t == tend["food kit",t]
            for j in DP
                updated[j] = false
            end
        end

        open(filename, "w") do f
        # ----------------------------------------------  demand ---------------------------------------------- #
            totalincrease = 0
            for j in DP
                if !updated[j] && rand() < 0.3 # j can only be updated once per food kit period
                    r = rand(-10:20)/100 # how much
                    c = "food kit"
                    ti = div(t-1,threshold[c]+1)*(threshold[c]+1)+1
                    howmuch = ceil(Int, abs(r)*sum(demand[c,j,t2] for t2 in ti:t))

                    if howmuch > 1
                        updateDemand(t, j, r, howmuch)
                        totalincrease += r > 0 ? howmuch : -howmuch
                        r = r > 0 ? 1 : -1
                        write(f, "$j\$$r\$$howmuch\n")
                        # println(t, " ", j, " ", r, " ", sum(demand[c,j,t2] for t2 in ti:t), " ", howmuch)
                    end
                end
            end

        # ----------------------------------------------  supply ---------------------------------------------- #
            if totalincrease > 1 && rand() < 0.5
                r = rand(30:100)/100 # how much
                totalincrease = ceil(Int, totalincrease*r)
                updateSupply(t, totalincrease)
                write(f, "$totalincrease\n")
            else
                write(f, "0\n")
            end
            # println(t, " ", totalincrease)
        end
    end

    # ----------------------------------------------  pi ---------------------------------------------- #
    filename = string(filepath, instance, "pied", t, ".txt")
    if !readPie(filename)
        up = false
        for i in N
            for j in N
                r = rand()
                if r < 0.1
                    pie[i,j,t] = 1
                    up = true
                elseif r > 0.9
                    pie[i,j,t] = 0
                    up = true
                end
            end
        end
        if rand() < 0.01
            up = true
            for i in N
                for j in N
                    if pie[i,j,t] == 0 || pie[i,j,t] > 1
                        for t2 in t:nT
                            if pie[i,j,t2] != 0
                                pie[i,j,t2] = 1.5
                            end
                        end
                    end
                end
            end
        end
        if up
            distance = floydWarshall()
            open(filename, "w") do f
                for (key,value) in pie
                    write(f, "$(key[1])\$$(key[2])\$$(key[3])\$$value\n")
                end
            end
        end
    end
end
