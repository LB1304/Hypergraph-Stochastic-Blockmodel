## Dataset description


```r
library(igraph)
bip_net <- read.csv("sandi_graph.csv", header = T, sep = " ")
```

We first construct the co-authorship network with the <tt>R</tt> package <tt>igraph</tt>
```r
Authors <- sort(unique(bip_net$Author))
Papers <- sort(unique(bip_net$Paper))
g <- graph.empty()
g <- add.vertices(g, nv = length(Authors), attr = list(name = Authors, type = rep(TRUE, length(Authors))))
g <- add.vertices(g, nv = length(Papers), attr = list(name = Papers, type = rep(FALSE, length(Papers))))
edgeListVec <- as.vector(t(as.matrix(data.frame(Author = bip_net$Author, Paper = bip_net$Paper))))
g <- add.edges(g, edgeListVec)
```

We inspect the structure of the the resulting graph in terms of connected components; we obtain 129 connected components, 
```r
cl <- components(g)
cl
```
We obtain 129 ; the size of the ... most populated is ... in the table below.





## Analysis with <tt>HyperSBM</tt>



## Estimation of the model
