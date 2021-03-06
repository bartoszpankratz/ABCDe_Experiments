#julia main.jl params_file

#params_file = ARGS[1]

using Base.Iterators
using Random
using ABCDGraphGenerator
using BenchmarkTools

function approx_mu(coms, d_min, d_max, γ, ξ)
    E_deg = ABCDGraphGenerator.get_ev(γ, d_min, d_max)
    W = E_deg * sum(coms)
    μ₀ = 1
    for c in coms
        μ₀ -= (c * E_deg/W)^2
    end
    ξ*μ₀
end

#Fixed Parameters:
seed = 33
max_iter = 1
isCL = false        
res_dir = "benchmark_res"

#Parameters:
n = 10 .^(range(4,9,step = 1))
ξs = range(0.0, 1.0, step = 0.05)
βs = [1.1, 1.5, 1.9]
γs = [2.1, 2.5, 2.9]
min_deg = [5,10]
is_local = [true, false]
threads = [1, 2, 4, 8, 16, 32];

mkpath(res_dir);
mkpath("logs")

Random.seed!(seed)

#params = sort(split.(readlines(params_file)),lt = (x,y)-> isless(x[1],y[1]));
params = sort(reshape(collect(product(n,ξs,βs,γs,min_deg)),:))

for (n, ξ, β, γ, min_deg) in params
    #get variables:
    n = parse(Int,n)
    ξ = parse(Float64,ξ)
    β = parse(Float64,β)
    γ = parse(Float64,γ)
    min_deg = parse(Int,min_deg)
    c_min = Int(0.005*n)               
    c_max = Int(0.2*n) 
    max_deg = Int(ceil(sqrt(n)))
    E_deg = ABCDGraphGenerator.get_ev(γ, min_deg, max_deg)
    
    (n >= 10^8 && min_deg == 10) && continue
    
    #sample degree sequence:
    deg_filename = "degs/deg_$(n)_$(min_deg)_$(max_deg)_$(γ).dat"
    if !isfile(deg_filename)
        degs = ABCDGraphGenerator.sample_degrees(γ, min_deg, max_deg, n, max_iter)
        open(io -> foreach(d -> println(io, d), degs), deg_filename, "w")
    end
    
    #sample communities:
    com_filename = "coms/com_$(n)_$(c_min)_$(c_max)_$(β).dat"
    if !isfile(com_filename)
        coms = ABCDGraphGenerator.sample_communities(β, c_min, c_max, n, max_iter)
        open(io -> foreach(d -> println(io, d), coms), com_filename, "w")
    else
        coms = parse.(Int,readlines(com_filename))
    end

    #approximate μ:
    μ = approx_mu(coms, min_deg, max_deg, γ, ξ)
    
    
    for islocal in is_local
        #ABCD algorithm:
        run(`julia get_ABCD.jl 
            $n $γ $min_deg $E_deg $deg_filename $β $com_filename 
            $μ $ξ $isCL $islocal $res_dir $seed`)
        @info "ABCD: $n $γ $min_deg $μ $ξ $islocal done!"
    end

    for t in threads
        #ABCDe algorithm:
        for islocal in is_local
            run(`julia --threads $t get_ABCDe.jl 
                $n $γ $min_deg $E_deg $deg_filename $β $com_filename 
                $μ $ξ $isCL $islocal $res_dir $seed`)
            @info "ABCDe: $t $n $γ $min_deg $μ $ξ $islocal done!"
        end

        #LFR algorithm:
        n <= 10^8 && run(`python get_LFR.py $t $n $γ $min_deg $E_deg $deg_filename $β $com_filename 
            $μ $ξ $res_dir $seed`)
        @info "LFR: $t $n $γ $min_deg $μ $ξ done!"
    end
end
@info "Experiment done!"

