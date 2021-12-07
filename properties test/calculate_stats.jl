# calculate_stats.jl algo n γ min_deg β ξ islocal res_dir no_iters

algo = ARGS[1]
n = parse(Int,ARGS[2])
γ = parse(Float64, ARGS[3])
min_deg = parse(Int, ARGS[4])
β = parse(Float64, ARGS[5])
ξ = parse(Float64, ARGS[6])
islocal = parse(Bool, ARGS[7])
res_dir = ARGS[8]
iters = parse(Int, ARGS[9])

using Graphs, StatsBase, Logging

degree_correlation(g) = [sum(sum(degree(g,nei) for nei in neighbors(g,n))/length(neighbors(g,n))
    for n in findall(x-> x == deg,degree(g)))/sum(degree(g).==deg)
        for deg in sort(unique(degree(g))) if deg != 0]

function correlation_exponent(g)
    x = log.(filter(x  -> x != 0, sort(unique(degree(g)))))
    y = log.(degree_correlation(g))
    # Fit the regression
    ([ones(length(x)) x] \ y)[2]
end

function centralities_stats(g, centralities, algo, flush)
    for cent in centralities
        name = "stats/$(algo)_$(Symbol(cent))_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).dat"
        open(name, "a") do io
            c = cent(g)
            avg_c, std_c = mean_and_std(c)
            Q1,Q2,Q3 = quantile(c,0.25), quantile(c,0.5), quantile(c,0.75)
            minc,maxc = minimum(c), maximum(c)
            res =[avg_c, std_c, Q1,Q2,Q3, minc, maxc]
            println(io, join(res,";"))
            if flush
            	cname = "stats/$(Symbol(cent))/$(algo)_$(Symbol(cent))_raw_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).dat"
            	cdf_func = ecdf(c)
            	open(cname, "w") do io
               	 for x in range(minimum(c), maximum(c), length = 100)
                    		println(io, x,"\t", cdf_func(x))
                	 end
            	end
            end
        end
    end
end

function shortest_paths(g, algo, flush)
    name = "stats/$(algo)_paths_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).dat"
    open(name, "a") do io
        #avg shortest path
        s = Int(ceil(0.1*nv(g))) 
        paths = vcat([dijkstra_shortest_paths(g, n).dists for n in rand(1:nv(g),s)]...)
        avg_shortest_path, std_shortest_path = mean_and_std(paths)
        Q1,Q2,Q3 = quantile(paths,0.25), quantile(paths,0.5), quantile(paths,0.75)
        minpath, maxpath = minimum(paths), maximum(paths)
        res =[avg_shortest_path, std_shortest_path, Q1,Q2,Q3, minpath, maxpath]
        println(io, join(res,";"))
        #flush all results to file
        if flush
        paths_name = "stats/paths/$(algo)_paths_raw_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).dat"
      	   cdf_func = ecdf(paths)
           open(paths_name, "w") do io
            	for path in sort(unique(paths))
               	 println(io, path, "\t", cdf_func(path))
            	end
           end
        end
    end
end

function participation_stats(g, partition, algo)
    name = "stats/$(algo)_participation_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).dat"
    open(name, "a") do io
        #Within module degree; hubs
        deg_int = [sum(partition[neighbors(g,n)] .== partition[n]) for n = 1:nv(g)]
        μA = [sum(deg_int[partition .== A])/sum(partition .== A) for A in sort(unique(partition))]
        σA = [sqrt(sum((deg_int[partition .== A] .- μA[A]).^2)/(sum(partition .== A))) for A in sort(unique(partition))]
        z = [(deg_int[n] - μA[partition[n]]) / σA[partition[n]] for n = 1:nv(g)] 
        hubs = sum(z .> 2.5)/nv(g)
        #participation coefficient
        part_coeffs = [degree(g,n) != 0 ? 1 - sum((sum(partition[neighbors(g,n)] .== A)/length(neighbors(g,n)))^2 
                for A in unique(partition)) : 0.0 for n = 1:nv(g)]
        avg_part_coeff = sum(part_coeffs)/nv(g)
        ultra_peripheral_nodes = sum(part_coeffs .< 0.05)/nv(g)
        peripheral_nodes = sum(0.05 .≤ part_coeffs .< 0.62)/nv(g)
        connector_nodes = sum(0.62 .≤ part_coeffs .< 0.8)/nv(g)
        kinless_nodes = sum(part_coeffs .≥ 0.8)/nv(g)
        provincial_hubs = sum(part_coeffs[z .> 2.5] .< 0.3) / (hubs*nv(g))
        connector_hubs = sum(0.3 .≤ part_coeffs[z .> 2.5] .< 0.75) / (hubs*nv(g))
        kinless_hubs = sum(part_coeffs[z .> 2.5] .≥ 0.75) / (hubs*nv(g))
        res =[avg_part_coeff, hubs, provincial_hubs, connector_hubs, kinless_hubs, 
            ultra_peripheral_nodes, peripheral_nodes, connector_nodes, kinless_nodes]
        println(io, join(res,";"))
    end
end

cores = [1,2,3]
centralities = [betweenness_centrality, closeness_centrality, degree_centrality,  pagerank]
(algo == "LFR" || algo == "LFRog") && (islocal = "")
name = "stats/$(algo)_stats_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).dat"
logname = "logs/main_log.dat"
mkpath("stats/knn")
mkpath("stats/paths")
for c in centralities
	mkpath("stats/$(Symbol(c))")
end

open(name, "w") do io
     itr = rand(1:iters)
     for iter in 1:iters
         try
         if algo != "LFRog"
                edge_list  = "$(res_dir)/$(algo)_$(iter)_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).dat"
                communities = "$(res_dir)/$(algo)_clusters_$(iter)_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).dat"
                start = 1
            else
                edge_list  = "$(res_dir)/$(algo)_$(iter)_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).nse"
                communities = "$(res_dir)/$(algo)_$(iter)_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).nmc"
                statsfile = "$(res_dir)/$(algo)_$(iter)_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).nst"
                start = 2
            end
            g = SimpleGraphFromIterator([Graphs.SimpleEdge(parse.(Int,e[1]),parse.(Int,e[2])) 
                for e in split.(readlines(edge_list))[start:end]])
            partition = [parse.(Int,e2) for (e1,e2) in split.(readlines(communities))]
            #modularity
            modularity = Graphs.modularity(g,partition)
            #clustering coefficients
            global_coeff = Graphs.global_clustering_coefficient(g)
            avg_local_coeff = sum(Graphs.local_clustering_coefficient(g,vertices(g)))/nv(g)
            res = [iter, modularity,global_coeff, avg_local_coeff]
            #number of nodes in k cores
            for k in cores
                push!(res,length(Graphs.k_core(g,min_deg + k))/nv(g))
            end
            #giant component
            giant_component = connected_components(g)[1]
            size = length(giant_component)/nv(g)
            volume = sum(degree(g,giant_component))/sum(degree(g))
            push!(res,size,volume)
            #robustness
            for p in [0.05,0.1,0.15]
                gcopy = deepcopy(g)
                rem_vertices!(gcopy, collect(1:ceil(Int, p * nv(g))))
                #size of graph
                push!(res, length(connected_components(gcopy)[1])/nv(g))
                #relative
                push!(res, length(connected_components(gcopy)[1])/nv(gcopy))
            end
            #internal edges
            comm_size = [sum(partition .== i) for i = 1:length(unique(partition))]
            internal_edges = zeros(length(unique(partition)))
            for v in 1:nv(g)
                sum = 0.0
                for nei in neighbors(g,v)
                    partition[v] == partition[nei] && (sum += 1)
                end
                internal_edges[partition[v]] += sum/degree(g,v)
            end
            push!(res,mean_and_std(internal_edges ./comm_size)...)
            #degree assortativity 
            correlation_coeff = assortativity(g)
            corr_exp = correlation_exponent(g)
            push!(res,correlation_coeff,corr_exp)
            println(io, join(res,";"))
            if iter == itr
            	knn_name = "stats/knn/$(algo)_knn_$(n)_$(ξ)_$(β)_$(γ)_$(min_deg)_$(islocal).dat"
            	knn = degree_correlation(g)
            	open(knn_name, "w") do io
               	 for (i,d) in enumerate([ d for d in sort(unique(degree(g))) if d != 0])
               	     println(io, d, "\t", knn[i])
               	 end
            	end
            end
            #participation stats
            participation_stats(g, partition, algo)
            shortest_paths(g, algo, iter == itr)
            centralities_stats(g, centralities, algo, iter == itr)
            #remove files
            rm(edge_list, force = true)
            rm(communities, force = true)
            algo == "LFRog" && rm(statsfile, force = true)
         catch err
	      logger = SimpleLogger(open(logname, "a"))
   		with_logger(logger) do
        		@warn exception = (err, stacktrace()), "  for params: $(algo) $(n) $(ξ) $(β) $(γ) $(min_deg) $(islocal)"
        		@warn ""
    		end
	  end
     
     end
end

    
        
