<h1 align="center">Analysis of co-authorship dataset</h1>
<p align="center"> <span style="font-size: 14px;"><em><strong>Luca Brusa &middot; Catherine Matias</strong></em></span> </p>
<br>


<h2>Dataset description</h2>

We analyze the Sandi co-authorship dataset, available at [this link](http://vlado.fmf.uni-lj.si/pub/networks/data/2mode/Sandi/Sandi.htm). The original dataset was extracted from the bibliography of the book "Product Graphs: Structure and Recognition" by _Imrich W._ and _Klav&#382;ar S._, and is given as a bipartite author/paper graph. The file sandi.net has been manually edited to create sandi_graph.csv that contains the list of edges in the bipartite graph Author/Paper.
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
g <- add.edges(g, edgeListVec)                                                                                  # Add the esges

A <- as_incidence_matrix(g)             # Obtain the incidence matrix
A <- A[!duplicated(A), ]                # Remove repeated papers
dim(A)
table(rowSums(A))/dim(A)[1]
```
The original dataset contains 274 papers and 314 authors.

We then remove all papers with more than 4 authors and all the ones with less than 2 authors. Moreover, we remove all authors that remain with no papers. Finally we remove repeated papers. 
```r
A <- A[-which(rowSums(A) > 4),]         # Remove papers with more than 4 authors 
A <- A[which(rowSums(A) > 1), ]         # Select only papers with more than 1 author
A <- A[, -which(colSums(A) == 0)]       # Remove corresponding columns authors with no papers
A <- A[!duplicated(A), ]                # Remove repeated papers
dim(A)
table(rowSums(A))/dim(A)[1]
```
The modified dataset has now 154 papers and 237 authors (70% of papers has 2 authors, 26% has 3 authors, and 4% has 4 authors).

We inspect the structure of the the resulting graph in terms of connected components: we obtain 59 connected components. We selected the main, containing 76 papers and 79 authors.
```r
g <- graph_from_incidence_matrix(A)
cl <- components(g)                                         # Create the connected components of the graph
vert_ids <- V(g)[cl$membership == which.max(cl$csize)]      # Select the main connected component
g_main <- induced_subgraph(g, vert_ids)                     # Build the corresponding bipartite graph

A_main <- as_incidence_matrix(g_main)                       # Obtain the incidence matrix
save(A_main, file = "./A_1cc.RData")
dim(A_main)
table(rowSums(A_main))/dim(A_main)[1]
```

We create the hypergraph in which nodes are authors and each hyperedge links the authors of a same paper. We save the hyperedges into the file HG_coauth_1cc.txt.
```r
hyperedges <- vector(mode = "list", length = nrow(A_main))
for (i in 1:nrow(A_main)) {
  ind <- unname(which(A_main[i, ] == 1))
  hyperedges[[i]] <- colnames(A_main)[ind]         # Store the hyperedges in a list
}

sink("./HG_coauth_1cc.txt")
for (i in 1:length(hyperedges)) {
  cat(paste(hyperedges[[i]], collapse = ","))      # Print the hyperedges to .txt file
  cat("\n")
}
sink()
```

---

<h2>Estimation with <tt>HyperSBM</tt></h2>

We fit the HSBM model on the resulting hypergraph using our package <tt>HyperSBM</tt>. We consider a number of latent groups ranging from 2 to 4, using two different initialization strategies: random and with soft spectral clustering;  ICL criterion selects $Q=2$ groups, and random initialization provides (here) the best results. (Further documentation for the <tt>HyperSBM</tt> package is available at [this link](https://github.com/LB1304/HyperSBM)).
```r
require(HyperSBM)
HG <- HyperSBM::import_hypergraph(file_name = "./HG_coauth_1cc.txt", method = "full")

set.seed(231)
for (q in 2:5) {
  res_rand <- HyperSBM::HSBM(Hypergraph = HG, Q = q, start = 0, model = 0, tol = 1e-5, maxit_VEM = 100, maxit_FP = 100, n_threads = 30, print = TRUE)
  res_ssc <- HyperSBM::HSBM(Hypergraph = HG, Q = q, start = 2, model = 0, tol = 1e-5, maxit_VEM = 100, maxit_FP = 100, n_threads = 30, print = TRUE)
  save(res_rand, res_ssc, file = paste0("res_coauth_1cc_Q", q, ".RData"))
}
```

---

<h2>Integrated Classification Likelihood and ELBO function plots</h2>

We compute the ICL criterion for each value of Q (Q = 2, 3, 4, 5) in order to select the optimal number of latent group. Th eplot shows that the highest ICL value is the one obtained with Q=2; therefore, relying on this criterion, we select 2 groups.
```r
load("./coauth_1cc_Q2_5_2init/res_coauth_1cc_Q2.RData")
if (res_rand$J > res_ssc$J){res_Q2 <- res_rand} else {res_Q2 <- res_ssc}
load("./coauth_1cc_Q2_5_2init/res_coauth_1cc_Q3.RData")
if (res_rand$J > res_ssc$J){res_Q3 <- res_rand} else {res_Q3 <- res_ssc}
load("./coauth_1cc_Q2_5_2init/res_coauth_1cc_Q4.RData")
if (res_rand$J > res_ssc$J){res_Q4 <- res_rand} else {res_Q4 <- res_ssc}
load("./coauth_1cc_Q2_5_2init/res_coauth_1cc_Q5.RData")
if (res_rand$J > res_ssc$J){res_Q5 <- res_rand} else {res_Q5 <- res_ssc}

plot(2:5, c(res_Q2$J, res_Q3$J, res_Q4$J, res_Q5$J), xlab = "Q", ylab = "ELBO")
plot(2:5, c(res_Q2$ICL, res_Q3$ICL, res_Q4$ICL, res_Q5$ICL), xlab = "Q", ylab = "ICL")
```

---

<h2>Analysis with <tt>HyperSBM<tt></h2>

Following ICL criterion, we select the model with $Q=2$ latent groups. We obtain a small group with only 8 authors, and a big group with the remaining 71 authors. 
```r
table(res_Q2$Z)                   # Frequency table in the 2 groups

small <- unname(which.min(table(res_Q2$Z)))
grpsmall_names <- names(colSums(A_main[, which(res_Q2$Z == small)]))  # Id of the authors in small group 1
grpsmall_names
grpsmall_ind <- which(colnames(A_main) %in% grpsmall_names)           # Column index in matrix A_main of the authors in small group 1
grpsmall_ind
```

We inspect the value of the variational parameter tau.
```r
round(res_Q2$tau, 2)
ind_ambiguous <- which(((0.2 < res_Q2$tau[, 1]) & (res_Q2$tau[, 1] < 0.8))) 
res_Q2$tau[ind_ambiguous, ]
```

We inspect the value of the variational parameter tau in the small group.
```r
grpsmall_tau <- res_Q2$tau[grpsmall_ind, ]
rownames(grpsmall_tau) <- grpsmall_names
round(grpsmall_tau, 2)
```

We inspect the value of the variational parameter tau in the big group.
```r
grpbig_tau <- res_Q2$tau[-grpsmall_ind, ]
rownames(grpbig_tau) <- setdiff(colnames(A_main), grpsmall_names)
round(grpbig_tau, 2)
```

We now analyze the characteristics of the small group. It results that among the 8 authors in small group 1:
* 6 of them have the highest number of distinct co-authors;
* 5 of them have the highest number of published papers.
```r
num_co_auth <- function(colonna) {          # Function to compute the number of co-authors of a given author (specified by its column index in matrix A_main)
  pap <- which(A_main[, colonna] != 0)
  co_auth <- c()
  for (p in pap) {
    co_auth <- c(co_auth, which(A_main[p, ] != 0))
  }
  return(length(unique(co_auth))-1)
}

table(unlist(lapply(1:ncol(A_main), num_co_auth)))   # Distribution of the number of distinct co-authors per author 
unlist(lapply(grpsmall_ind, num_co_auth))            # Number of distinct co-authors for the 8 authors in small group 1
unlist(lapply(ind_ambiguous, num_co_auth))           # Same for ambiguous authors

table(colSums(A_main))                               # Degree distribution of authors
colSums(A_main[, which(res_Q2$Z == small)])          # Degree of the 8 authors in small group 1
colSums(A_main[, ind_ambiguous])                     # Same for ambiguous authors
```

We also inspect the values of the estimated probabilities of hyperedges occurrance. The highest probabilities are obtained, in general, for probabilities of hyperedges between nodes from different groups. This shows that neither tha first nor the second group are communities.
```r
res_Q2$B          # Estimates for the probability of occurrence of an hyperedge
```

---

<h2>Comparison with Hypergraph Spectral Clustering</h2>

We first look at the eigenvalues of the Laplacian function; the largest gap method selects 15 groups.
```r
L <- HG$Laplacian
lambda <- eigen(L)$values
sum(lambda < 10^(-6))               # We recover the number of cc by looking at the zero eigenvalues
plot(lambda)
which.max(abs(lambda[length(lambda):2]-lambda[(length(lambda)-1):1]))       # Selected number of groups by largest spectral gap method  
```

Then we compare our approach with the spectral clustering algorithm proposed in Ghoshdastidar and Dukkipati (2017), using Q=2 groups.
```r
SpectralClust <- function (L, Q, n) {              # Estimation function for Spectral Clustering
  X <- as.matrix(eigen(L)$vectors[, (n - Q + 1):n])
  X <- X/apply(X,1,norm,type="2")
  km <- kmeans(x = X, centers = Q, nstart = 100)
  
  return(km$cluster)
}

res_sc <- SpectralClust(L = HG$Laplacian, Q = 2, n = HG$Num_nodes)

table(res_sc)               # Frequency table in the groups
small <- unname(which.min(table(res_sc)))

grpsmall_sc_names <- colnames(A_main)[which(res_sc == small)]       # Id of the authors in small group 1
grpsmall_sc_ind <- which(colnames(A_main) %in% grpsmall_sc_names)   # Column index in matrix A_main of the authors in small group 1

grp2_sc_names <- colnames(A_main)[which(res_sc != small)]    # Id of the authors in large group 2
grp2_sc_ind <- which(colnames(A_main) %in% grp2_sc_names)    # Column index in matrix A_main of the authors in large group 2

hist(unlist(lapply(grpsmall_sc_ind, num_co_auth)))           # Number of distinct co-authors for the authors in the smaller group
table(unlist(lapply(grpsmall_sc_ind, num_co_auth)))
hist(unlist(lapply(grp2_sc_ind, num_co_auth)))         # Number of distinct co-authors for the authors in the larger group
table(unlist(lapply(grp2_sc_ind, num_co_auth)))

hist(colSums(A_main[, which(res_sc == small)]))        # Degree of the authors in the first group
table(colSums(A_main[, which(res_sc == small)]))
hist(colSums(A_main[, which(res_sc != small)]))        # Degree of the authors in the large group
table(colSums(A_main[, which(res_sc != small)]))
```

We also compute the B probabilities corresponding to spectral clustering groups.
```r
M <- HG$Max_size
n <- HG$Num_nodes
all_Latents <- vector(mode = "list", length = M-1)
for (m in 2:M) {
  all_Latents[[m - 1]] <- RcppAlgos::comboGeneral(v = 2, m = m, repetition = TRUE)
}
list_H.Edges <- HG$List_of_H.Edges
Y <- compute_Y(n = n, M = M, list_edges = list_H.Edges)

tau_sc <- matrix(0, ncol = 2, nrow = length(res_sc))
for (i in 1:nrow(tau_sc)) {
  tau_sc[i, res_sc[i]] <- 1
}

B_sc <- HyperSBM::compute_B(n = n, M = M, tau = tau_sc, Y = Y, all_latents = all_Latents, model = 0, n_threads = 1)
B_sc
```

---

<h2>Comparison with BipartiteSBM</h2>

We then analyze the bipartite graph of authors/papers with the <tt>R</tt> package <tt>sbm</tt> with option `bipartite`. 
```r
require(sbm)
res_bp <- estimateBipartiteSBM(A_main)       # Estimation function for Bipartite SBM
table(res_bp$memberships$col)                # Frequency table in the groups
# 2 groups of authors are selected
small <- unname(which.min(table(res_bp$memberships$col)))

# Look at posterior probabilities in small group
round(res_bp$probMemberships$col[,small],2)
# ambiguous nodes 
bp_ind_amb <- which(((0.2<res_bp$probMemberships$col[,small]) &(res_bp$probMemberships$col[,small]<0.8))) 
bp_ind_amb

grpsmall_bp_names <- colnames(A_main)[which(res_bp$memberships$col == small)]    # Id of the authors in small group 1
grpsmall_bp_ind <- which(colnames(A_main) %in% grpsmall_bp_names)                # Column index in matrix A_main of the authors in small group 1
grpsmall_bp_ind
grpsmall_bp_names

table(unlist(lapply(1:ncol(A_main), num_co_auth)))   # Distribution of the number of distinct co-authors per author 
unlist(lapply(grpsmall_bp_ind, num_co_auth))         # Number of distinct co-authors for the authors in the smaller group
unlist(lapply(bp_ind_amb, num_co_auth))              # Same for ambiguous authors

table(colSums(A_main))                                      # Degree distribution of authors
colSums(A_main[, which(res_bp$memberships$col == small)])   # Degree of the authors in the first group
colSums(A_main[, bp_ind_amb])                               # Degree of the authors in the ambiguous group

grpsmall_bp_tau <- res_bp$probMemberships$col[grpsmall_ind, ]
rownames(grpsmall_bp_tau) <- grpsmall_names
round(grpsmall_bp_tau, 2)         # The first 4 rows are the ones in small group according to BipartiteSBM

grp2_bp_tau <- res_bp$probMemberships$col[-grpsmall_ind, ]
rownames(grp2_bp_tau) <- setdiff(colnames(A_main), grpsmall_names)
round(grp2_bp_tau, 2)

# Connectivity parameters in bipartite SBM
res_bp$connectParam

# Compute B values for Bipartite clusters 
gp_bp <- apply(res_bp$probMemberships$col,1,which.max)
tau_bp <- matrix(0, ncol = 2, nrow = length(gp_bp))
for (i in 1:nrow(tau_bp)) {
  tau_bp[i, gp_bp[i]] <- 1
}

B_bp <- HyperSBM::compute_B(n = n, M = M, tau = tau_bp, Y = Y, all_latents = all_Latents, model = 0, n_threads = 1)
B_bp
```




