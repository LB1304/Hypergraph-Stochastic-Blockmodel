##########
# Julia code to use Kaminski's method on Line clustering data
#########

# Uncomment those lines if you require pakages installation
#=
using Pkg
Pkg.add("SimpleHypergraphs")
Pkg.add("DelimitedFiles")
Pkg.add("Clustering")
=#

using SimpleHypergraphs
using Clustering

####### KAMINSKI AUXILIARY ########

function ari(x,y)
    evaluations = randindex(x, y)
    ari = evaluations[1]
    return ari
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

# Function to read data and construct Hypergraph object
function do_file(name::String)
    f = open(name)
    line = readline(f) # here we read the first line of the file - should be 0,1,2,...,n-1
    h = SimpleHypergraphs.Hypergraph{Bool,Int}(0,0) # create empty hypergraph

    for v_meta in parse.(Int,(split(line,",")))
        # print(v_meta)
        add_vertex!(h,v_meta=v_meta) # add each vertex from first line in the hypergraph
    end
    
    for line in eachline(f) # then continue reading the lines, but starting from the second one !
        x = parse.(Int,(split(line,",")))
        inds = x #.+ 1  # each hyperedge is encoded by numbering nodes from 1 to n
        # print(inds)
        add_hyperedge!(h;vertices=Dict(inds .=> true))
    end
    close(f)
    return h
end


function find_comms(h, nreps::Int=1000; ha = SimpleHypergraphs.HypergraphAggs(h))
    # Initialize modularity
    best_modularity = 0
    # Initialize partition with all vertices in its own part
    comms = [Set(i) for i in 1:nhv(h)]  #nhv return the number of vertices in the hypergraph h
    mod_history = Vector{Float64}(undef, nreps)
    
    for rep in 1:nreps
        
        # Randomly select an hyperedge (using simplified stochastic algorithm version)
        he = rand(1:nhe(h))
        vers = collect(keys(getvertices(h, he))) #Get vertices returns vertices from a hypergraph h for a given hyperedge he
        c = deepcopy(comms)   #copy mutuable object
        
        # Compute partition obtained when merging all parts in 'c' touched by 'he'
        i0 = find_first(c, vers)
        max_i = length(c)
        i_cur = i0
        #print(i_cur)
        #print(max_i)
        while i_cur < max_i
            i_cur += 1
            if length(intersect(c[i_cur],vers)) > 0
                union!(c[i0], c[i_cur])
                c[i_cur]=c[max_i]
                max_i += -1
            end
        end
        resize!(c,max_i)
        
        # Compute modularity
        m = SimpleHypergraphs.modularity(h, c, ha)
        if m > best_modularity
            best_modularity = m
            
            # update partition
            comms = c
        end
                
        mod_history[rep] = best_modularity
    end
    return (bm=best_modularity, bp=comms, mod_history=mod_history)
end

function find_first(c::Array{Set{Int}}, vals)
    for i in 1:length(c)
        for v in vals
            v in c[i] && return i
        end
    end
end


# Transform partition Z outputed by Kaminski into a vector of length number of nodes
function transform_partition(Z, n)
    
    Z_hat= zeros(Int64, n)
    # associate nodes to clusters: singletons stay in different communities
    for i in range(1,length(Z))
        for j in Z[i]
            Z_hat[j]=i
        end
    end

    k=0
    # add singleton clusters to nodes that are not in some partitions
    for i in range(1,length(Z_hat))
        if Z_hat[i]==0
            k += 1
            Z_hat[i]=length(Z)+k
        end
    end
    return Z_hat
end

##################
## end of auxiliary functions
###################

REP=99
Ari = zeros(Float64, REP+1)
Qhat = zeros(Float64, REP+1)


for rep in 0:REP

    # Switch here between 2lines and 3lines cases
    #dataname = "data_2linecluster/rep_" * string(rep) * "/"
    #n=100
    dataname = "data_3linecluster/rep_" * string(rep) * "/"
    n=150
    
    # Read isolated nodes and remove them
    isolated = read_isolated(dataname)

    # Read true clusters
    Z_true = read_clusters(dataname)

    # Create Hyperedges_add.txt file adding a line describing set of nodes
    firstline = join(string.(1:n), ",") # we keep all nodes here
    hyperedges = read("$dataname/hyperedges.txt", String)
    concat = firstline * "\n" * hyperedges
    write("$dataname/hyperedges_add.txt", concat )
  
    # Construct hypergraph and apply algorithm
    h = do_file("$dataname/hyperedges_add.txt")
    # Set the number of iterations to 3*tot nb hyperedges
    niter = size(h)[2]*3    #  nb of iterations over hyperedges
    ha = SimpleHypergraphs.HypergraphAggs(h)
    # modularity is unstable, we compute it 20 times and keep best    
    result_mod = -1e-12
    for k in 1:20
        res = find_comms(h,niter;ha)
        if res[1] > result_mod
            global Z = res[2]
            result_mod = res[1]
        end
    end  
    #
    Z_hat = transform_partition(Z,n)
    # Remove isolated nodes
    deleteat!(Z_hat, sort(isolated))
    Qhat[rep+1] = length(unique(Z_hat))
    
    # performance measures
    Ari[rep+1]=ari(Z_hat,Z_true)

    # End of repetition
    println("repetition ", rep, " finished")
    
end



#### Save results ####

# Switch here between 2lines and 3lines cases
#=
open("Kaminski_ARI_2lines.txt", "w") do f
    for ari in Ari
        write(f, "$ari\n")
    end
end

open("Kaminski_Qhat_2lines.txt", "w") do f
    for q in Qhat
        write(f, "$q\n")
    end
end
=#

#######
# Switch here between 2lines and 3lines cases

open("Kaminski_ARI_3lines.txt", "w") do f
    for ari in Ari
        write(f, "$ari\n")
    end
end

open("Kaminski_Qhat_3lines.txt", "w") do f
    for q in Qhat
        write(f, "$q\n")
    end
end
