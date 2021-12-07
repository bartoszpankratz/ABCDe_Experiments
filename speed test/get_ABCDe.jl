# gen_ABCDe.jl n γ min_deg E_deg degs  β coms μ ξ isCl islocal path seed

max_iter = 1000
n = parse(Int,ARGS[1])
γ = parse(Float64, ARGS[2])
min_deg = parse(Int, ARGS[3])
E_deg = parse(Float64, ARGS[4])
degs = ARGS[5]
β = parse(Float64, ARGS[6])
coms = ARGS[7]
μ = parse(Float64, ARGS[8])
ξ = parse(Float64, ARGS[9])
isCL = parse(Bool, ARGS[10])
islocal = parse(Bool, ARGS[11])
mu = μ
xi = ξ
islocal == false ? mu = nothing : xi = nothing
res_dir = ARGS[12]

ABCDename = "$(res_dir)/ABCDe_$(Threads.nthreads())_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).dat"
logname = "logs/ABCDe_$(Threads.nthreads())_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).log"

(isfile(ABCDename) || isfile(logname)) && exit()

using Logging
using BenchmarkTools
using Random
using ABCDeGraphGenerator

Random.seed!(parse(Int, ARGS[13]))

@eval BenchmarkTools macro btimed(args...)
           _, params = BenchmarkTools.prunekwargs(args...)
           bench, trial, result = gensym(), gensym(), gensym()
           trialmin, trialallocs = gensym(), gensym()
           tune_phase = hasevals(params) ? :() : :($BenchmarkTools.tune!($bench))
           return esc(quote
               local $bench = $BenchmarkTools.@benchmarkable $(args...)
               $BenchmarkTools.warmup($bench)
               $tune_phase
               local $trial, $result = $BenchmarkTools.run_result($bench)
               local $trialmin = $BenchmarkTools.minimum($trial)
               $result, $BenchmarkTools.time($trialmin)/1e9
           end)
       end

try 
    n_edges, time = BenchmarkTools.@btimed begin
        #_ = ABCDeGraphGenerator.sample_degrees($γ, $min_deg, $Int(ceil(sqrt(n))), $n, $max_iter)
        #_ = ABCDeGraphGenerator.sample_communities($β, $Int(0.005*n), $Int(0.2*n), $n, $max_iter)
        p = ABCDeGraphGenerator.ABCDParams($parse.(Int,readlines(degs)), $parse.(Int,readlines(coms)), 
                                            $mu, $xi, $isCL, $islocal)
        edges, _ = ABCDeGraphGenerator.gen_graph(p)
        length(edges)
        end
    @assert n_edges == sum(parse.(Int,readlines(degs)))/2 
    res = [n, ξ, μ, β, γ, min_deg, E_deg, n_edges, islocal, Threads.nthreads(), time]
    open(io -> foreach(r -> println(io, r), res), ABCDename, "w")
catch err
    logger = SimpleLogger(open(logname, "w+"))
    with_logger(logger) do
        @warn exception = (err, stacktrace()), "  for params: $(n) $(ξ) $(β) $(γ) $(min_deg) $(islocal)"
    end
end
