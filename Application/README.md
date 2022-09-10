<h1 align="center">Analysis of co-authorship dataset</h1>
<br>
<em align="center">Luca Brusa &middot; Catherine Matias</em>

<br>


<h2>Dataset description</h2>

We analze a co-authorship dataset available at [this link](http://vlado.fmf.uni-lj.si/pub/networks/data/2mode/Sandi/Sandi.htm). The original dataset was extracted from the bibliography of the book "Product Graphs: Structure and Recognition" by _Imrich W._ and _Klav&#382;ar S._, and is given as a bipartite author/paper graph.
```r
require(igraph)
bip_net <- read.csv("sandi_graph.csv", header = T, sep = " ")
```

We construct the co-authorship bipartite network with the <tt>R</tt> package <tt>igraph</tt>. 
```r
Authors <- sort(unique(bip_net$Author))     # Id of authors
Papers <- sort(unique(bip_net$Paper))       # Id of papers
g <- graph.empty()
g <- add.vertices(g, nv = length(Authors), attr = list(name = Authors, type = rep(TRUE, length(Authors))))      # Add the first layer of nodes (authors)
g <- add.vertices(g, nv = length(Papers), attr = list(name = Papers, type = rep(FALSE, length(Papers))))        # Add the second layer of nodes (papers)
edgeListVec <- as.vector(t(as.matrix(data.frame(Author = bip_net$Author, Paper = bip_net$Paper))))
g <- add.edges(g, edgeListVec)                                                                                  # Add the edges
```

We inspect the structure of the the resulting graph in terms of connected components; we obtain 129 connected components, the main containing 86 authors and 167 papers.
```r
cl <- components(g)                                         # Create the connected components of the graph
vert_ids <- V(g)[cl$membership == which.max(cl$csize)]      # Select the main connected component
g_main <- induced_subgraph(g, vert_ids)                     # Build the corresponding bipartite graph
```

We create the hypergraph in which nodes are authors and each hyperedge links the authors of a same paper. We choose to discard the 62 papers published by a unique author and the unique paper published by 5 co-authors (hyperedges of sizes 1 and 5). This results into a hypergraph with 83 authors and 104 hyperedges of sizes between 2 and 4 (71%, 27%, and 2% of sizes 2, 3, and 4 respectively).
```r
A <- as_incidence_matrix(g_main)                    # Obtain the incidence matrix
A <- A[-1, ]                                        # Remove the paper with 5 authors
A_mod <- A[, -which(colSums(A) == 0)]
A_mod <- A_mod[which(rowSums(A_mod) > 1), ]         # Select only papers with more than 1 author

hyperedges <- vector(mode = "list", length = nrow(A_mod))
for (i in 1:nrow(A_mod)) {
    ind <- unname(which(A_mod[i, ] == 1))
    hyperedges[[i]] <- colnames(A_mod)[ind]         # Store the hyperedges in a list
}
hyperedges <- unique(hyperedges)                    # Remove repeated hyperedges

sink("./HG_coauth.txt")
for (i in 1:length(hyperedges)) {
    cat(paste(hyperedges[[i]], collapse = ","))     # Print the hyperedges to .txt file
    cat("\n")
}
sink()
```

---

<h2>Analysis with <tt>HyperSBM</tt></h2>

We fit the HSBM model on the resulting hypergraph using our package <tt>HyperSBM</tt>. We consider a number of latent groups ranging from 2 to 4, using two different initialization strategies: random and with soft spectral clustering;  ICL criterion selects $Q=2$ groups, and random initialization provides (here) the best results. (Further documentation for the <tt>HyperSBM</tt> package is available at [this link](https://github.com/LB1304/HyperSBM)).
```r
require(HyperSBM)
HG <- HyperSBM::import_hypergraph(file_name = "./HG_coauth.txt", method = "full")

for (q in 2:4) {
  res_rand <- HyperSBM::hSBM_par(Hypergraph = HG, Q = q, start = 0, model = 0, tol = 1e-5, maxit_VEM = 100, maxit_FP = 100, n_threads = 30, print = TRUE)
  res_ssc <- HyperSBM::hSBM_par(Hypergraph = HG, Q = q, start = 2, model = 0, tol = 1e-5, maxit_VEM = 100, maxit_FP = 100, n_threads = 30, print = TRUE)
  save(res_rand, res_ssc, file = paste0("res_coauth_Q", q, ".RData"))
}
```

Following ICL criterion, we select the model with $Q=2$ latent groups. We obtain a small group with only 6 authors, and a big group with the remaining 77 authors. 
```r
load("./res_coauth_Q2.RData")
table(res_rand$Z)                   # Frequency table in the 2 groups
small <- unname(which.min(table(res_rand$Z)))
```

The 6 authors in small group 1 have the highest number of distinct co-authors.
```r
num_co_auth <- function(colonna) {              # Function to compute the number of co-authors of a given author (specified by its column index in matrix A_mod)
  pap <- which(A_mod[, colonna] != 0)
  co_auth <- c()
  for (p in pap) {
    co_auth <- c(co_auth, which(A_mod[p, ] != 0))
  }
  return(length(unique(co_auth))-1)
}

grp1_names <- names(colSums(A_mod[, which(res_rand$Z == small)]))       # Id of the authors in small group 1
grp1_ind <- which(colnames(A_mod) %in% grp1_names)                      # Column index in matrix A_mod of the authors in small group 1

table(unlist(lapply(1:ncol(A_mod), num_co_auth)))   # Distribution of the number of distinct co-authors per author 
unlist(lapply(grp1_ind, num_co_auth))               # Number of distinct co-authors for the 6 athors in small group 1
```

In this small group we have 4 of the 5 authors that wrote the highest number of papers (highest hyperedge degree); this group also contains an author with smaller degree.
```r
table(colSums(A_mod))                               # Degree distribution of authors
colSums(A_mod[, which(res_rand$Z == small)])        # Degree of the 6 authors in small group 1
```

We also inspect the values of the estimated probabilities of hyperedges occurrance. The highest probabilities are obtained, in general, for probabilities of hyperedges between nodes from different groups. This shows that neither tha first nor the second group are communities.
```r
res_rand$B          # Estimates for the probability of occurrance of an hyperedge
```

---

<h2>Comparison with two other methods</h2>

We first compare our approach with the spectral clustering algorithm proposed in Ghoshdastidar and Dukkipati (2017).
```r
SpectralClust <- function (L, Q, n) {                       # Estimation function for Spctral Clustering
    X <- as.matrix(eigen(L)$vectors[, (n - Q + 1):n])
    X <- X/rowSums(X)
    km <- kmeans(x = X, centers = Q, nstart = 100)
    
    return(km$cluster)
}
res_sc <- SpectralClust(L = HG$Laplacian, Q = 2, n = HG$Num_nodes)

table(res_sc)               # Frequency table in tha groups
small <- unname(which.min(table(res_sc)))

grp1_sc_names <- colnames(A_mod)[which(res_sc == small)]    # Id of the authors in small group 1
grp1_sc_ind <- which(colnames(A_mod) %in% grp1_sc_names)    # Column index in matrix A_mod of the authors in small group 1

colSums(A_mod[, which(res_sc == small)])        # Degree of the authors in the first group
unlist(lapply(grp1_sc_ind, num_co_auth))        # Number of distinct co-authors for the authors in the smaller group
```

We then analyze the bipartite graph of authors/papers with the <tt>R</tt> package <tt>sbm</tt> with option `bipartite`. 
```r
require(sbm)
res_bp <- estimateBipartiteSBM(A_mod)       # Estimation function for Bipartite SBM
table(res_bp$memberships$col)               # Frequency table in tha groups
small <- unname(which.min(table(res_bp$memberships$col)))

grp1_bp_names <- colnames(A_mod)[which(res_bp$memberships$col == small)]    # Id of the authors in small group 1
grp1_bp_ind <- which(colnames(A_mod) %in% grp1_bp_names)                    # Column index in matrix A_mod of the authors in small group 1

colSums(A_mod[, which(res_bp$memberships$col == small)])    # Degree of the authors in the first group
unlist(lapply(grp1_bp_ind, num_co_auth))                    # Number of distinct co-authors for the authors in the smaller group
```




