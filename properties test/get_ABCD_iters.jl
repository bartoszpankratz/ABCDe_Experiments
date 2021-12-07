# gen_ABCDe.jl n γ min_deg degs β coms μ ξ isCl islocal res_dir no_iters

n = parse(Int,ARGS[1])
γ = parse(Float64, ARGS[2])
min_deg = parse(Int, ARGS[3])
degs = ARGS[4]
β = parse(Float64, ARGS[5])
coms = ARGS[6]
μ = parse(Float64, ARGS[7])
ξ = parse(Float64, ARGS[8])
isCL = parse(Bool, ARGS[9])
islocal = parse(Bool, ARGS[10])
mu = μ
xi = ξ
islocal == false ? mu = nothing : xi = nothing
res_dir = ARGS[11]
iters = parse(Int, ARGS[12])

using BenchmarkTools
using Logging
using Random
using ABCDGraphGenerator

for iter = 1:iters
    Random.seed!(iter + 21)
    ABCDname = "$(res_dir)/ABCD_$(iter)_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).dat"
    comname = "$(res_dir)/ABCD_clusters_$(iter)_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).dat"
    logname = "logs/ABCD_$(iter)_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).log"
    try 
        p = ABCDGraphGenerator.ABCDParams(parse.(Int,readlines(degs)), parse.(Int,readlines(coms)), 
                                            mu, xi, isCL, islocal)
        edges, clusters = ABCDGraphGenerator.gen_graph(p)

        open(ABCDname, "w") do io
            for (a, b) in sort!(collect(edges))
                println(io, a, "\t", b)
            end
        end

        open(comname, "w") do io
            for (i, c) in enumerate(clusters)
                println(io, i, "\t", c)
            end
        end
    catch err
        logger = SimpleLogger(open(logname, "w+"))
        with_logger(logger) do
            @warn @warn exception = (err, stacktrace()), "  for params: $(n) $(μ) $(β) $(γ) $(min_deg) $(islocal)"
        end
    end
end