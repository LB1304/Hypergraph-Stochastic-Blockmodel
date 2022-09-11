<h1 align="center">Simulation studies</h1>
<p align="center"> <span style="font-size: 14px;"><em><strong>Luca Brusa &middot; Catherine Matias</strong></em></span> </p>
<br>

<h2>Parameters and groups estimation</h2>

We simulate hypergraphs from the HSBM model considering a simplified latent structure. We consider $Q=2$ latent groups with priori probabilities equal to 0.4 and 0.6 respectively. The largest size $M$ of hyperedges is set to 3. Four different values are examined for the number of nodes: $n=50, 100, 150, 200$.
```r
n_vec <- c(50, 100, 150, 200)
for (i in 1:length(n_vec)) {
  for (n_rep in 1:10) {
    HyperSBM::sample_Hypergraph(n = n_vec[i], M = 3, Q = 2, pi = c(0.4, 0.6), alpha = 0.25, beta = 0.35, file_name = paste0("HG", n_rep))
    
    HG <- HyperSBM::import_Hypergraph(file_name = paste0("./HG", n_rep, ".txt"), method = "full")
    
    out_full <- HyperSBM::HSBM(Hypergraph = HG, Q = 2, start = 2, model = 0, tol = 1e-6, maxit_VEM = 25, maxit_FP = 25, n_threads = 6)
    out_aff  <- HyperSBM::HSBM(Hypergraph = HG, Q = 2, start = 2, model = 1, tol = 1e-6, maxit_VEM = 25, maxit_FP = 25, n_threads = 6)
    
    save(out_full, out_aff, file = paste0("./HG", n_rep, "_results.RData"))
  }
}
```



<h2>Model selection</h2>
