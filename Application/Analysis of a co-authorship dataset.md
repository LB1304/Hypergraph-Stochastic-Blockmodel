## Dataset description 

We analze a co-authorship dataset available at [this link](http://vlado.fmf.uni-lj.si/pub/networks/data/2mode/Sandi/Sandi.htm). The original dataset was extracted from the bibliography of the book "Product Graphs: Structure and Recognition" by _Imrich W._ and _Klav&#382;ar S._, and is given as a bipartite author/paper graph.
```r
library(igraph)
bip_net <- read.csv("sandi_graph.csv", header = T, sep = " ")
```

We construct the co-authorship network with the <tt>R</tt> package <tt>igraph</tt>. 
```r
Authors <- sort(unique(bip_net$Author))
Papers <- sort(unique(bip_net$Paper))
g <- graph.empty()
g <- add.vertices(g, nv = length(Authors), attr = list(name = Authors, type = rep(TRUE, length(Authors))))
g <- add.vertices(g, nv = length(Papers), attr = list(name = Papers, type = rep(FALSE, length(Papers))))
edgeListVec <- as.vector(t(as.matrix(data.frame(Author = bip_net$Author, Paper = bip_net$Paper))))
g <- add.edges(g, edgeListVec)
```

We inspect the structure of the the resulting graph in terms of connected components; we obtain 129 connected components, the main containing 86 authors and 167 papers.
```r
cl <- components(g)
vert_ids <- V(g)[cl$membership == which.max(cl$csize)]
g_main <- induced_subgraph(g, vert_ids)
```

We create the hypergraph in which nodes are authors and each hyperedge links the authors of a same paper. We choose to discard the 62 papers published by a unique author and the unique paper published by 5 co-authors (hyperedges of sizes 1 and 5). This results into a hypergraph with 83 authors and 104 hyperedges of sizes between 2 and 4 (71%, 27%, and 2% of sizes 2, 3, and 4 respectively).
```r
A <- as_incidence_matrix(g_main)
A <- A[-1, ]
A_mod <- A[, -which(colSums(A) == 0)]
A_mod <- A_mod[which(rowSums(A_mod) > 1), ]

hyperedges <- vector(mode = "list", length = nrow(A_mod))
for (i in 1:nrow(A_mod)) {
    ind <- unname(which(A_mod[i, ] == 1))
    hyperedges[[i]] <- colnames(A_mod)[ind]
}
hyperedges <- unique(hyperedges)

sink("./HG_coauth.txt")
for (i in 1:length(hyperedges)) {
    cat(paste(hyperedges[[i]], collapse = ","))
    cat("\n")
}
sink()
```


## Analysis with <tt>HyperSBM</tt>

We fit the HSBM model on the resulting hypergraph using our package <tt>HyperSBM</tt>. We consider 
```r

```


## Estimation of the model
