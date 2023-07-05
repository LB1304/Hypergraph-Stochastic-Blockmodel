<h1 align="center">Simulation studies</h1>
<p align="center"> <span style="font-size: 14px;"><em><strong>Luca Brusa &middot; Catherine Matias</strong></em></span> </p>
<br>

<h2>Parameters and groups estimation</h2>

We simulate hypergraphs from the HSBM model considering a simplified latent structure. We explore three different scenarios:
- Communities: high intra-groups and low inter-groups connection probabilities (&alpha; > &beta;);
- Disassortative: low intra-groups and high inter-groups connection probabilities (&alpha; < &beta;);
- Erdös-Rényi-like : very similar intra-groups and inter-groups connection probabilities (&alpha; $\approxeq$ &beta;).

For each scenario values of $\alpha^{(2)}$ and $\beta^{(2)}$ decrease with increasing $n$; specific value are summarized in the following Table. Moreover, $\alpha^{(3)} = \alpha^{(2)} / n$ and $\beta^{(3)} = \beta^{(2)} / n$.

<table>
  <tr>
    <td></td><td colspan="2">Scenario A</td><td colspan="2">Scenario B</td><td colspan="2">Scenario C</td>
  </tr>
  <tr>
    <td>n</td><td>a2</td><td>b2</td><td>a2</td><td>b2</td><td>a2</td><td>b2</td>
  </tr>
</table>


We consider $Q=2$ latent groups with priori probabilities equal to 0.4 and 0.6 respectively. The largest size $M$ of hyperedges is set to 3. Four different values are examined for the number of nodes: $n=50, 100, 150, 200$.
```r
n_vec <- c(50, 100, 150, 200)
for (i in 1:length(n_vec)) {
  for (n_rep in 1:10) {
    # Draw samples from the HSBM with "Communities scenario" (for the other scenarios it is enough to modify the values of alpha and beta)
    HyperSBM::sample_Hypergraph(n = n_vec[i], M = 3, Q = 2, pi = c(0.4, 0.6), alpha = 0.7, beta = 0.3, file_name = paste0("HG", n_rep))
    
    # Import the hypergraph
    HG <- HyperSBM::import_Hypergraph(file_name = paste0("./HG", n_rep, ".txt"), method = "full")
    
    # Estimate the HSBM model
    out_full <- HyperSBM::HSBM(Hypergraph = HG, Q = 2, start = 2, model = 0, tol = 1e-6, maxit_VEM = 25, maxit_FP = 25, n_threads = 8)
    
    save(out_full, out_aff, file = paste0("./HG", n_rep, "_results.RData"))
  }
}
```



<h2>Model selection</h2>

We simulate hypergraphs from the HSBM with $Q=3$ assuming the simplified latent structure (in the "Communities") scenario. 
The largest size $M$ of hyperedges is set to 3. Two different values are examined for the number of nodes: $n=100, 200$.
In the estimation phase, we consider a number of latent groups ranging from 1 to 5.
```r
for (n_rep in 1:50) {
  # Draw samples
  HyperSBM::sample_Hypergraph(n = 100, M = 3, Q = 3, pi = rep(1/3, 3), alpha = 0.7, beta = 0.3, file_name = paste0("HG", n_rep))
}

for (n_rep in 1:50) {
  # Import the hypergraph
  HG <- HyperSBM::import_Hypergraph(file_name = paste0("./HG", n_rep, ".txt"))
  
  # Estimate the HSBM modle for Q ranging from 1 to 5
  res <- vector(mode = "list", length = 5)
  for (q in 1:5) {
    res[[q]] <- HyperSBM::HSBM(Hypergraph = HG, Q = q, start = 2, model = 0, tol = 1e-6, maxit_VEM = 25, maxit_FP = 25, n_threads = 8)
  }

  save(res, file = paste0("./HG", n_rep, "_results.RData"))
}
```




