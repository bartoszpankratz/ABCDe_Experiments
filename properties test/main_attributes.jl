#julia main.jl 

using Base.Iterators
using Random
using ABCDGraphGenerator
using StatsBase

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
seed = 21
max_iter = 1
isCL = false        
res_dir = "graphs"
no_iters = 30

#Parameters:
n = 10^4
ξs = range(0.0, 1.0, step = 0.05)
βs = [1.1, 1.5, 1.9]
γs = [2.1, 2.5, 2.9]
#min_deg = [1,5,10]
min_deg = 5
is_local = [true, false]
threads = [1, 4, 16, 32];

mkpath(res_dir);
mkpath("logs")
mkpath("stats")

Random.seed!(seed)

params = sort(reshape(collect(product(n,ξs,βs,γs,min_deg)),:))

for (n, ξ, β, γ, min_deg) in params
    c_min = Int(ceil(0.005*n))               
    c_max = Int(ceil(0.2*n)) 
    max_deg = Int(ceil(sqrt(n)))
    #sample degree sequence:
    deg_filename = "degs/deg_$(n)_$(min_deg)_$(max_deg)_$(γ).dat"
    if !isfile(deg_filename)
        degs = ABCDGraphGenerator.sample_degrees(γ, min_deg, max_deg, n, 1)
        open(io -> foreach(d -> println(io, d), degs), deg_filename, "w")
    end
    
    #sample communities:
    com_filename = "coms/com_$(n)_$(c_min)_$(c_max)_$(β).dat"
    if !isfile(com_filename)
        coms = ABCDGraphGenerator.sample_communities(β, c_min, c_max, n, 1)
        open(io -> foreach(d -> println(io, d), coms), com_filename, "w")
    else
        coms = parse.(Int,readlines(com_filename))
    end

    #approximate μ:
    μ = approx_mu(coms, min_deg, max_deg, γ, ξ)
    E_deg = ABCDGraphGenerator.get_ev(γ, min_deg, max_deg)
    
    #LFR algorithm:
    run(`python get_LFR_iters.py $n $γ $min_deg $deg_filename $β $com_filename 
          $μ $ξ $res_dir $no_iters`)
    run(`julia calculate_stats.jl LFR $n $γ $min_deg  $β $ξ false $res_dir $no_iters`)
    @info "LFR: $n $γ $min_deg $μ $ξ done!"
    
    #LFR OG algorithm:
    for iter = 1:no_iters
        try
            run(`timeout 5 ./LFR-Benchmark_UndirWeightOvp/lfrbench_udwov -N $n -k $(Int(ceil(E_deg)))  -maxk $max_deg -muw $μ
                 -t1 $β -t2 $γ -minc $c_min -maxc $c_max -dname $deg_filename -cname $com_filename
                -name $(res_dir)/LFRog_$(iter)_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_`)
        catch
            nothing
        end
    end
    run(`julia calculate_stats.jl LFRog $n $γ $min_deg  $β $ξ false $res_dir $no_iters`)
    @info "LFR OG: $n $γ $min_deg $μ $ξ done!"
    
    for islocal in is_local
        #ABCD algorithm:
        run(`julia get_ABCD_iters.jl 
            $n $γ $min_deg $deg_filename $β $com_filename 
            $μ $ξ $isCL $islocal $res_dir $no_iters`)
        run(`julia calculate_stats.jl ABCD $n $γ $min_deg  $β $ξ $islocal $res_dir $no_iters`)
        @info "ABCD: $n $γ $min_deg $μ $ξ $islocal done!"
    end
end
@info "Experiment done!"
