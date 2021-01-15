using JuMP, CPLEX, DataStructures, DelimitedFiles, Statistics, Printf#, MathOptFormat

const EXACT = 1
const CG = 2
const PS = 3
const KNOWFUTURE = false

const BYPRIORITY = 1
const BYDEMAND = 2

const B = 1000000
const ZERO = 0.0001
const NTHREADS = 8
const SLEEP = 10
const PERCENTAGE = 10
const PERIODS_PS = 3

instance = ARGS[1]
# const METHOD = parse(Int, ARGS[3])
# const FIRSTCHOICE = parse(Int, ARGS[4])
# const PSEQ = parse(Bool, ARGS[5])
# const nV = parse(Float64, ARGS[6])
# cities = readdlm("/home/pcandida/punim0040/instances/cidadesSP.txt")
# auxDistance = readdlm("/home/pcandida/punim0040/instances/distanciasSP.txt")/1000
# filepath = "/home/pcandida/punim0040/instances/deterministic"
#instance = "151"
const METHOD = PS
const FIRSTCHOICE = BYPRIORITY
const PSEQ = true
const nV = parse(Int, instance) > 5 ? 2.5 : 3.5
cities = readdlm("/Users/pcandida/Dropbox/unimelb/planilhas/brasil/cidadesSP.txt")
auxDistance = readdlm("/Users/pcandida/Dropbox/unimelb/planilhas/brasil/distanciasSP.txt")/1000
filepath = "/Users/pcandida/Dropbox/unimelb/instancesBRd/"

println(filepath, " ", instance, " ", METHOD, " ", FIRSTCHOICE, " ", PSEQ, " ", nV)
const DEBUG = true #instance == "67"

const o = "BAURU, SP"#["APIAI, SP", "BAURU, SP", "PRESIDENTE PRUDENTE, SP", "SAO PAULO, SP", "REGISTRO, SP", "RIBEIRAO PRETO, SP", "TAUBATE, SP"]# Array{String}(readdlm(string(filepath, instance, "FN.txt"))[:,1]) # candidates
populationDP = readdlm(string(filepath, instance, "DP.txt"))
DP = Array{String}(populationDP[:,1])
populationDP = Array{Integer}(populationDP[:,2])

const N = sort([o; DP])
const No = ["theta"; N]
nDP = length(DP)
n = nDP + 1

const IMPROVEMENT = 0.98
const IMPROVEMENTB = 0.70
const MAXARCS = nDP < 26 ? 36 : 18
# const TL = nDP < 26 ? 20 : 40
const TIME_LIMIT = nDP * (nDP+1) / 2  #testar tb com 30s

distance = Dict{Any,Float64}()
for i in 1:645
    for j in 1:645
        if cities[i] in N && cities[j] in N
            distance[cities[i],cities[j]] = auxDistance[i,j]
        end
    end
end
max_distance = maximum(auxDistance)

const nC = 3
const L = ["relief itens", "food"]
const C = ["food kit", "cleaning kit", "mattress"]
# const compatibility = Dict{Any,Int8}()
# compatibility["food kit", "relief itens"] = 1
# compatibility["cleaning kit", "relief itens"] = 1
# compatibility["mattress", "relief itens"] = 1
# compatibility["food kit", "food"] = 1
# compatibility["cleaning kit", "food"] = 1
# compatibility["mattress", "food"] = 1
const priority = Dict{String,Integer}("food kit" => 15, "cleaning kit" => 2, "mattress" => 3)
const weight = Dict{String,Float64}("food kit" => 5, "cleaning kit" => 1, "mattress" => 3)
# const volume = Dict{String,Float64}("food kit" => 0.001, "cleaning kit" => 0.001, "mattress" => 0.1)

const nT = 15
const T = [1:nT;]
const working_hours = 12
const threshold = Dict{String,UInt8}()
threshold["food kit"] = 4 # threshold["food kit"] = period => order expire for periodic demand
threshold["cleaning kit"] = 9
threshold["mattress"] = 14
const Tend = Dict{String,Array{Integer}}()
auxset = Set()
for c in C
    Tend[c] = Array{Integer}(undef,0)
    tmax = threshold[c]+1
    while tmax < nT
        push!(Tend[c], tmax)
        push!(auxset, tmax)
        tmax += threshold[c]+1
    end
    push!(Tend[c], nT)
end
push!(auxset, nT)
const Tendall = sort(collect(auxset))
const tend = Dict{Any, Integer}()
for t in T
    for c in C
        tend[c,t] = min(div(t+threshold[c],threshold[c]+1)*(threshold[c]+1),nT)
    end
end

# supply - http://www.spherehandbook.org/en/water-supply-standard-1-access-and-water-quantity/
const peopleDP = Dict{String,Integer}()
for j in 1:nDP
    peopleDP[DP[j]] = trunc(Integer, populationDP[j]*PERCENTAGE/100)
end
const total_demand = sum(values(peopleDP))
# println("% of people affected = ", PERCENTAGE, ", number of people affected = ", total_demand)

const supply = Dict{Any,Float64}()
const demand = Dict{Any,Float64}()
const urgentdemand = Dict{Any,Float64}()
for t in T
    supply["food kit",t]     = (t-1) % 5 == 0 ? trunc(Integer, total_demand) : 0
    supply["cleaning kit",t] = t == 1 || t == 11 ? trunc(Integer, total_demand/4) : 0
    supply["mattress",t]     = t == 1 ? trunc(Integer, total_demand) : 0

    for (j,pop) in peopleDP
        urgentdemand["food kit",j,t]     = (t-1) % 5 == 0 ? pop : 0
        urgentdemand["cleaning kit",j,t] = t == 1 || t == 11 ? trunc(Integer, pop/4) : 0
        urgentdemand["mattress",j,t]     = t == 1 ? pop : 0
        demand["food kit",j,t]     = urgentdemand["food kit",j,t]
        demand["cleaning kit",j,t] = urgentdemand["cleaning kit",j,t]
        demand["mattress",j,t]     = urgentdemand["mattress",j,t]
    end
end

# capacity FN
# http://www.logcluster.org/blog/logistics-cluster-bolsters-ongoing-humanitarian-response-fiji-mobile-storage-units
# MSU is 24x10x5 meters, or a total of 240 m², with a capacity to store 840 cubic meters of food and Non Food Items (NFIs).
# 840 m3 = 70% of 1200 m3
const capacityF = 840
# opening costs - not used

# H-36 Caracal
# https://www.airbus.com/helicopters/military-helicopters/heavy/h225m.html Range = 452 NM = 837 Km, 2551 L * 2.07 = R$ 5280 => 6.31
# https://www.airbushelicoptersinc.com/products/H225-specifications.asp
# http://tecnodefesa.com.br/fab-recebe-helicoptero-h-36-caracal-na-versao-operacional/
# http://estudio.folha.uol.com.br/brasil-que-voa/2017/05/1886629-preco-do-combustivel-de-aviacao-no-brasil-e-46-maior-do-que-nos-eua.shtml
# Can we assume re-fuelling will take part with loading?

# https://estradao.estadao.com.br/caminhoes/os-10-caminhoes-mais-vendidos-em-2017/
# Volkswagen Delivery 8.160 Mercedes-Benz Sprinter 415 CDI
# https://www.noticiasautomotivas.com.br/sprinter/
# 15-16 mpg $1 diesel/L
const K = Array{String}(undef,0)
const nvehicles = Dict{Any,Integer}()
const vehicle = Dict{Any,Bool}()
const vtypes = ["MB", "VW"]
const vmodel = Dict{String,String}()
const vload = Dict{String,Integer}("MB" => 2040, "VW" => 8150)#, "helicopter" => 5406
const vloadtime = Dict{String,Float16}("MB" => 0.7, "VW" => 1.5, "helicopter" => 0.5)
const vcost = Dict{String,Integer}("MB" => 100, "VW" => 280, "helicopter" => 1000) # R$ per day
const vfuel = Dict{String,Float64}("MB" => 2.1016*0.12, "VW" => 2.1016*0.15, "helicopter" => 6.31) # R$ per km
const vspeed = Dict{String,Float64}("MB" => 70, "VW" => 70) # km/h

totalWeight = 0
for c in C
    global totalWeight += weight[c] * sum(supply[c,t] for t in T)
end
const nK = totalWeight < ZERO ? 0 : totalWeight/(2*nT) <= 2040 ? 1 : 2*ceil(Int,totalWeight/(nV*(sum(values(vload))*nT-8150)))
aux = ceil(Int, nK/2)
for k in 1:nK
    push!(K, string("k",k))
    vmodel[K[k]] = totalWeight/(2*nT) <= 2040 || k > aux ? "MB" : "VW"
end

for t in T
    for k in vtypes
        nvehicles[k,t] = 0
    end
end
for t in T
    for k in K
        vehicle[k,t] = vmodel[k] != "VW" || t != 1
        nvehicles[vmodel[k],t] += vmodel[k] == "VW" && t == 1 ? 0 : 1
    end
end
# println("nvehicles = ", sum(nvehicles[vtype,2] for vtype in vtypes), " nK = ", nK)
#=aux = readdlm(string(filepath, instance, "vehicle.txt"), '$')
for i in 1:trunc(Int, length(aux)/3)
    vehicle[aux[i,1],aux[i,2]] = aux[i,3]
end
open(string(filepath, instance, "vehicle.txt"), "w") do f
    for (key,value) in vehicle
        write(f, "$(key[1])\$$(key[2])\$$value\n")
    end
end=#

pie = Dict{Any,Float32}()
filename = string(filepath, instance, "pied.txt")
if readPie(filename)
    local aux = readdlm(filename, '$')
    for i in 1:trunc(Int, length(aux)/4)
        pie[aux[i,1],aux[i,2],aux[i,3]] = aux[i,4]
    end
else
    pie = getPi()
    open(filename, "w") do f
        for (key,value) in pie
            write(f, "$(key[1])\$$(key[2])\$$(key[3])\$$value\n")
        end
    end
end
distance = floydWarshall()

nP = min(trunc(Int, working_hours/(minDistance(distance)*2/maximum(values(vspeed)) + minimum(values(vloadtime)))) + 2, 6)
P = collect(1:nP)
# @printf("T = %d, working_hours = %d, P = %d\n", nT, working_hours, nP)

const wfairness = 100

# const S = Array{String}(undef, 0)
# const Sk = Dict{String,Array{String}}("MB" => Array{String}(undef, 0), "VW" => Array{String}(undef, 0))
# const columncost = Dict{Any,Float64}()
# const columndeltaQ = Dict{Any, Float64}()
# const columxvalue = Dict{Any, Float64}()

# release memory and force garbage collector to run
populationDP = auxv = nothing
GC.gc()
