##########
# Julia code to use Chodrow's method on Line clustering data
#########

# Uncomment those lines if you require pakages installation
#=
using Pkg
Pkg.add("HyperModularity")
Pkg.add("SimpleHypergraphs")
Pkg.add("Clustering")
=#

using HyperModularity
using SimpleHypergraphs
using Clustering


##################
## Auxiliary functions from Chodrow et al. 
###################

function ari(x,y)
    evaluations = randindex(x, y)
    ari = evaluations[1]
    return ari
end

# Hypergraph has nodes encoded from 1 to n 
function read_hypergraph_edges(dataname::String, maxsize::Int64=25, minsize::Int64=2)
    E = Dict{Integer, Dict}()
    filename = join([dataname,"hyperedges.txt"])
    open(filename) do f
        for line in eachline(f)
            edge = [parse(Int64, v) for v in split(line, ',')]
            sort!(edge)
            if minsize <= length(edge) <= maxsize
                sz = length(edge)
                if !haskey(E, sz)
                    E[sz] = Dict{}()
                end
                E[sz][edge] = 1
            end
        end
    end
    return E
end

function read_clusters(dataname::String)
    labels = Int64[]
    filename = join([dataname,"true_clusters.txt"])
    open(filename) do f
        for line in eachline(f)
            push!(labels, parse(Int64, line))
        end
    end
    return labels
end


function read_isolated(dataname::String)
    isolated = Int64[]
    filename = join([dataname,"isolated.txt"])
    open(filename) do f
        for line in eachline(f)
            push!(isolated, parse(Int64, line))
        end
    end
    return isolated
end


##################
## end of auxiliary functions
###################


REP=99
AriAON = zeros(Float64, REP+1)
AriSymm = zeros(Float64, REP+1)
Qhat_AON = zeros(Float64, REP+1)
Qhat_Symm = zeros(Float64, REP+1)

for rep in 0:REP 
    
    # Switch here between 2lines and 3lines cases
    dataname = "data_2linecluster/rep_" * string(rep) * "/"
    n=100
    #dataname = "data_3linecluster/rep_" * string(rep) * "/"
    #n=150
    
    # Read isolated nodes
    isolated = read_isolated(dataname)

    # Read true clusters
    Z_true = read_clusters(dataname)
    
    N = 1:n
    E = read_hypergraph_edges(dataname)
    D = zeros(Int64, n)
    for (sz, edges) in E
        for (e, _) in edges
            D[e] .+= 1
        end
    end
    maxedges = maximum(keys(E))
    for k in 1:maxedges
        if !haskey(E, k)
            E[k] = Dict{}()
        end
    end
    H = hypergraph(N, E, D)
    #print(H)
    

    ############ AON HMLL ################
    #### Start with CLIQUELOUVAIN ####
    
    start = "cliquelouvain" 
    Z_AON = Simple_AON_Louvain(H, startclusters = start)
    
    # delete isolated nodes from result
    deleteat!(Z_AON, sort(isolated))

    Qhat_AON[rep+1] = length(unique(Z_AON))
    #print("Number of clusters found: $Q")
    AriAON[rep+1] = ari(Z_AON,Z_true)
    
    ############ SYMMETRIC HMLL ################

    kmax = maximum(keys(E))
    function ω(p, α)
        k = sum(p)
        return sum(p)/sum((p .* (1:length(p)).^α[k])) / n^(α[kmax+k]*k)
    end
    Ω = partitionIntensityFunction(ω, kmax)
    alpha = vcat(repeat([1.0], kmax), 0.1*(1:kmax)) # initial guess
    
    nrounds = 10
    for i ∈ 1:nrounds
        α = alpha
        global Z_Symm = HyperModularity.SuperNodeLouvain(H, Ω; α, verbose = false)
        alpha = HyperModularity.learnParameters(H, Z_Symm, Ω, α; tol = 1e-8)
        Q = round(Float64(HyperModularity.modularity(H, Z_Symm, Ω; α)), digits = 2)
    end
    
    # delete isolated nodes from result
    deleteat!(Z_Symm, sort(isolated))
    
    Qhat_Symm[rep+1] = length(unique(Z_Symm))
    #print("Number of clusters found: $Q")
    AriSymm[rep+1] = ari(Z_Symm , Z_true)

    # End of repetition
    println("repetition ", rep, " finished")
    
end  



#### Save results ####

# Switch here between 2lines and 3lines cases
open("Chodrow_AON_ARI_2lines.txt", "w") do f
    for ari in AriAON  
        write(f, "$ari\n")
    end
end

open("Chodrow_AON_Qhat_2lines.txt", "w") do f
    for q in Qhat_AON
        write(f, "$q\n")
    end
end  

open("Chodrow_Symm_ARI_2lines.txt", "w") do f
    for ari in AriSymm
        write(f, "$ari\n")
    end
end

open("Chodrow_Symm_Qhat_2lines.txt", "w") do f
    for q in Qhat_Symm
        write(f, "$q\n")
    end
end


#######
# Switch here between 2lines and 3lines cases

open("Chodrow_AON_ARI_3lines.txt", "w") do f
    for ari in AriAON  
        write(f, "$ari\n")
    end
end  

open("Chodrow_AON_Qhat_3lines.txt", "w") do f
    for q in Qhat_AON
        write(f, "$q\n")
    end
end  

open("Chodrow_Symm_ARI_3lines.txt", "w") do f
    for ari in AriSymm
        write(f, "$ari\n")
    end
end

open("Chodrow_Symm_Qhat_3lines.txt", "w") do f
    for q in Qhat_Symm
        write(f, "$q\n")
    end
end


