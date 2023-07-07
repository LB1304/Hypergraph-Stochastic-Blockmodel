<h1 align="center">Simulation studies under HSBM</h1>
<p align="center"> <span style="font-size: 14px;"><em><strong>Luca Brusa &middot; Catherine Matias</strong></em></span> </p>
<br>

<h2>Parameters and groups estimation</h2>

We simulate sparse hypergraphs from the HSBM considering a simplified latent structure. We explore three different scenarios:
- (A) Communities: high intra-groups and low inter-groups connection probabilities (&alpha; > &beta;);
- (B) Disassortative: low intra-groups and high inter-groups connection probabilities (&alpha; < &beta;);
- (C) Erdös-Rényi-like : very similar intra-groups and inter-groups connection probabilities (&alpha; $\approxeq$ &beta;).

For each scenario values of $\alpha^{(2)}$ and $\beta^{(2)}$ (for pairwise intra-group and inter-groups interactions, respectively) decrease with increasing $n$ (such that $n\alpha^{(2)}$ and $n\beta^{(2)}$ remain constant); specific values are summarized in the following Table. Moreover, $\alpha^{(3)}, \beta^{(3)}$ (for size-3 intra-group and inter-groups interactions, respectively) are such that $\alpha^{(3)} = \alpha^{(2)} / n$ and $\beta^{(3)} = \beta^{(2)} / n$.

<div align="center">
<table>
  <tr>
    <td></td><td colspan="2">Scenario A</td><td colspan="2">Scenario B</td><td colspan="2">Scenario C</td>
  </tr>
  <tr>
    <td>$n$</td><td>$\alpha^{(2)}$</td><td>$\beta^{(2)}$</td><td>$\alpha^{(2)}$</td><td>$\beta^{(2)}$</td><td>$\alpha^{(2)}$</td><td>$\beta^{(2)}$</td>
  </tr>

  <tr>
    <td>50</td> <td>0.7000</td> <td>0.3000</td> <td>0.3000</td> <td>0.7000</td> <td>0.2500</td> <td>0.3500</td>
  </tr>
  <tr>
    <td>100</td> <td>0.3500</td> <td>0.1500</td> <td>0.1500</td> <td>0.3500</td> <td>0.1250</td> <td>0.1720</td>
  </tr>
  <tr>
    <td>150</td> <td>0.2300</td> <td>0.1000</td> <td>0.1000</td> <td>0.2300</td> <td>0.0800</td> <td>0.1200</td>
  </tr>
  <tr>
    <td>200</td> <td>0.1750</td> <td>0.0750</td> <td>0.0750</td> <td>0.1750</td> <td>0.0625</td> <td>0.0875</td>
  </tr>
</table>
</div>


We consider $Q=2$ latent groups with probabilities equal to 0.4 and 0.6 respectively. The largest size $M$ of hyperedges is set to 3. Four different values are examined for the number of nodes: $n=50, 100, 150, 200$.
```r
# Draw a sample from the m-Aff HSB sub-model under Scenario A and considering 100 nodes
# Saves it in the "HG.txt" file
HyperSBM::sample_Hypergraph(n = 50, M = 3, Q = 2, pi = c(0.4, 0.6), alpha = 0.35, beta = 0.15, file_name = "HG")
    
# Import the hypergraph
HG <- HyperSBM::import_Hypergraph(file_name = "./HG.txt", method = "full")
    
# Estimate the full HSB model using the Soft Spectral Clustering initialization
out_full <- HyperSBM::HSBM(Hypergraph = HG, Q = 2, start = 2, model = 0, tol = 1e-6, maxit_VEM = 25, maxit_FP = 25, n_threads = 8)

# Estimate the m-Aff HSB sub-model using the Absolute Spectral Clustering initialization
out_aff_m <- HyperSBM::HSBM(Hypergraph = HG, Q = 2, start = 3, model = 2, tol = 1e-6, maxit_VEM = 25, maxit_FP = 25, n_threads = 8)  
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




