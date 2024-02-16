using JuMP, CPLEX, DataStructures, DelimitedFiles, Statistics, MathOptInterface, Combinatorics

const B = 100000000
const ZERO = 0.0001
const NTHREADS = 8
const SLEEP = 1

const CHONLY = 0
const RELAXFIX = 1
const FIXOPT = 2
const PS = 3
const EXACT = 4
const DC = 5
const GRASP = 6
# const PETAL = 7

const ITUPDATE = 100
alpha = 0.0

const BYPRIORITY = 1
const BYDEMAND = 2

instance = ARGS[1]
const ninstance = parse(Int, instance)
const METHOD = parse(Int, ARGS[2])
const RUN_LS = false #parse(Bool, ARGS[7])
const RERUN_LS = false # compared, makes no difference
const FIRSTCHOICE = parse(Int, ARGS[3])
const EQCONSTR = parse(Bool, ARGS[4])
const psupply = parse(Float64, ARGS[5]) 
const pvehicles = parse(Int, ARGS[6])
const multivehicle = true #parse(Bool, ARGS[7])
const multicommodity = true # parse(Bool, ARGS[8])
const roadpi = true # parse(Bool, ARGS[9])
const DEBUG = true #instance == "67"
const KNOWFUTURE = false #parse(Bool, ARGS[8])
const DYNAMIC_ANALYSIS = false #parse(Bool, ARGS[9])
const wfairness = parse(Int, ARGS[7])
const wtransp = parse(Int, ARGS[8])

cities = readdlm("/home/pcandida/punim0040/instances/cidadesSP.txt")
coordinates = readdlm("/home/pcandida/punim0040/instances/coordenadasSP.txt")
auxDistance = readdlm("/home/pcandida/punim0040/instances/distanciasSP.txt")/1000
filepath = "/home/pcandida/punim0040/instances/deterministic/"
# cities = readdlm("/Users/pcandida/Dropbox/unimelb/planilhas/brasil/cidadesSP.txt")
# coordinates = readdlm("/Users/pcandida/Dropbox/unimelb/planilhas/brasil/coordenadasSP.txt")
# auxDistance = readdlm("/Users/pcandida/Dropbox/unimelb/planilhas/brasil/distanciasSP.txt")/1000
# filepath = "/Users/pcandida/Dropbox/unimelb/instancesBRd/"

#["APIAI, SP", "BAURU, SP", "PRESIDENTE PRUDENTE, SP", "SAO PAULO, SP", "REGISTRO, SP", "RIBEIRAO PRETO, SP", "TAUBATE, SP"]# Array{String}(readdlm(string(filepath, instance, "FN.txt"))[:,1]) # candidates
const o = ninstance < 190 ? "BAURU, SP" : "SAO PAULO, SP"
populationDP = readdlm(string(filepath, instance, "DP.txt"))
DP = Array{String}(populationDP[:,1])
populationDP = Array{Integer}(populationDP[:,2])

const N = sort([o; DP])
const No = ["theta"; N]
nDP = length(DP)
n = nDP + 1

const IMPROVEMENT = 0.99
# const IMPROVEMENTB = 0.70
const TL = 4
MY_TIME_LIMIT = (div(nDP * (nDP+1), 3) + 30)# / (METHOD == GRASP ? TL : 3) #(ninstance > 15 && ninstance < 200) || (ninstance > 215) ? 4 : 2)
const WINDOW = 2
const PERCENTAGE = ninstance < 190 ? 10 : 3
println(filepath, " ", instance, " ", PERCENTAGE, "% supply: ", psupply, " vehicles: ", pvehicles, " CH: ", FIRSTCHOICE, " IH: ", METHOD, " eq: ", EQCONSTR, " future: ", KNOWFUTURE, " time limit: ", MY_TIME_LIMIT)

coordinate = Dict{Any,Tuple{Float64,Float64}}()
for i in 1:645
    if coordinates[i,1] in N
        coordinate[coordinates[i,1]] = (coordinates[i,2],coordinates[i,3])
    end
end

const nC = multicommodity ? 3 : 1
const L = ["relief itens", "food"]
const C = multicommodity ? ["food kit", "cleaning kit", "mattress"] : ["food kit"]
# const compatibility = Dict{Any,Int8}()
# compatibility["food kit", "relief itens"] = 1
# compatibility["cleaning kit", "relief itens"] = 1
# compatibility["mattress", "relief itens"] = 1
# compatibility["food kit", "food"] = 1
# compatibility["cleaning kit", "food"] = 1
# compatibility["mattress", "food"] = 1
const priority = Dict{String,Integer}("food kit" => 7, "cleaning kit" => 2, "mattress" => 3)
const priorityo = Dict{String,Integer}("food kit" => 15, "cleaning kit" => 2, "mattress" => 3)
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
const total_demand = psupply * PERCENTAGE/100 * sum(populationDP)
supply, demand = getSupplyDemand()
println("% of people affected = ", PERCENTAGE, ", number of people affected on first week = ", sum(demand[C[1],j,t] for j in DP for t in 1:5),
    ", supply = ", sum(supply[C[1],t] for t in 1:5))

# capacity FN
# http://www.logcluster.org/blog/logistics-cluster-bolsters-ongoing-humanitarian-response-fiji-mobile-storage-units
# MSU is 24x10x5 meters, or a total of 240 m², with a capacity to store 840 cubic meters of food and Non Food Items (NFIs).
# 840 m3 = 70% of 1200 m3
# const capacityF = 840
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
# caminhão toco https://blog.texaco.com.br/ursa/tipos-de-caminhoes-e-capacidades/
const unloadtime = 0.25
const K = Array{String}(undef,0)
const nvehicles = Dict{Any,Integer}()
const vehicle = Dict{Any,Bool}()
const vtypes = multivehicle ? ["MB", "VW"] : ["toco"] #"MB", "VW"]
const vmodel = Dict{String,String}()
const vload = Dict{String,Integer}(multivehicle ? ("MB" => 2000, "VW" => 8150) : "toco" => 7000)
const vloadtime = Dict{String,Float16}(multivehicle ? ("MB" => 0.25, "VW" => 0.6) : "toco" => 0.5) # per trip, not per stop, "helicopter" => 0.5
const vcost = Dict{String,Integer}(multivehicle ? ("MB" => 200, "VW" => 300) : "toco" => 200) # R$ per day , "helicopter" => 1000
const vfuel = Dict{String,Float64}(multivehicle ? ("MB" => 7/6, "VW" => 7/5) : "toco" => 7/4) # R$ per km , "helicopter" => 6.31
const vspeed = Dict{String,Float64}(multivehicle ? ("MB" => 90, "VW" => 90) : "toco" => 80) #, "VW" => 70) # km/h

#nKdict = Dict([(1,2), (2,2), (3,2), (4,1), (5,3), (6,3), (7,3), (8,2), (9,3), (10,1), (11,4), (12,2), (13,2), (14,2), (15,2), (16,5), (17,3), (18,3), (19,2), (20,3), (21,4), (22,3), (23,7), (24,7), (25,5), (26,8), (27,6), (28,5), (29,6), (30,4), (31,4), (32,8), (33,7), (34,6), (35,6), (191,1), (192,1), (193,1), (194,1), (195,1), (196,1), (197,1), (198,1), (199,1), (200,1), (201,2), (202,3), (203,2), (204,3), (205,3), (206,2), (207,3), (208,4), (209,3), (210,3), (211,3), (212,3), (213,2), (214,5), (215,4), (216,5), (217,4), (218,3), (219,4), (220,4), (221,4), (222,5), (223,7), (224,7), (225,4)])
nKdict = Dict([(1,1), (2,1), (3,2), (4,1), (5,1), (6,1), (7,1), (8,2), (9,1), (10,1), (11,4), (12,1), (13,1), (14,2), (15,1), (16,3), (17,2), (18,3), (19,2), (20,2), (21,4), (22,3), (23,5), (24,5), (25,3), (26,7), (27,5), (28,4), (29,4), (30,4), (31,3), (32,7), (33,6), (34,5), (35,6), (191,1), (192,1), (193,1), (194,1), (195,1), (196,1), (197,1), (198,1), (199,1), (200,1), (201,2), (202,2), (203,2), (204,3), (205,2), (206,2), (207,3), (208,3), (209,2), (210,3), (211,3), (212,3), (213,2), (214,4), (215,3), (216,3), (217,4), (218,3), (219,4), (220,4), (221,4), (222,4), (223,5), (224,5), (225,4)])
nK = get(nKdict, ninstance, 0) - pvehicles
#nK = max(round(Int,(ninstance > 190 ? 1 : 1.5)*18.5*total_demand/(2.5*sum(values(vload))*nT)), (multivehicle && psupply > 0.99 && ninstance > 200 ? 2 : 1)) - 
	#pvehicles + (multivehicle && (ninstance in [204,207,208,210,211,212] || ninstance > 213) ? 1 : 0)
# aux = ninstance > 211 ? ceil(Int, nK/2) : floor(Int, nK/2)
aux = floor(Int, nK/2)
for k in 1:nK
    push!(K, string("k",k))
    vmodel[K[k]] = multivehicle ? ((nK == 1 && ninstance in [194,195]) || k <= aux ? "VW" : "MB") : "toco"
end

for t in T
    for k in vtypes
        nvehicles[k,t] = 0
    end
end
for t in T
    for k in K
        vehicle[k,t] = true #vmodel[k] != "VW" || t != 1
        nvehicles[vmodel[k],t] += 1 # vmodel[k] == "VW" && t == 1 ? 0 : 1
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

distance = Dict{Any,Float64}()
if !roadpi
for i in 1:645
    for j in 1:645
        if cities[i] in N && cities[j] in N
            for t in T
                distance[cities[i],cities[j],t] = auxDistance[i,j]
            end
        end
    end
end
max_distance = maximum(auxDistance)
else
filename = string(filepath, instance, "pied.txt")
pie = Dict()
if !readPie(filename)
    println("ERROR: could not find ", filename)
    exit(1)
    # for i in 1:645
    #     for j in 1:645
    #         if cities[i] in N && cities[j] in N
    #             for t in T
    #                 distance[cities[i],cities[j]] = auxDistance[i,j]
    #             end
    #         end
    #     end
    # end
    # pie = getPi()
    # open(filename, "w") do f
    #     for (key,value) in pie
    #         write(f, "$(key[1])\$$(key[2])\$$(key[3])\$$value\n")
    #     end
    # end
    # distance = floydWarshall()
end
end

nP = trunc(Int, working_hours/(minimum(distance[o,j,t] for j in DP for t in T)*2/maximum(values(vspeed)) + minimum(values(vloadtime)) + unloadtime)) + 2
P = collect(1:nP)
println("multicommodity ", multicommodity, ", pi ", roadpi, ", multivehicle ", multivehicle, ", pvehicles ", pvehicles, ", nK = ", nK, ", nP = ", nP)


# const S = Array{String}(undef, 0)
# const Sk = Dict{String,Array{String}}("MB" => Array{String}(undef, 0), "VW" => Array{String}(undef, 0))
# const columncost = Dict{Any,Float64}()
# const columndeltaQ = Dict{Any, Float64}()
# const columxvalue = Dict{Any, Float64}()

# release memory and force garbage collector to run
populationDP = auxv = nKdict = nothing
GC.gc()
flush(stdout)
GC.gc()
