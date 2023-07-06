<h1 align="center">Line clustering through hypergraphs</h1>
<p align="center"> <span style="font-size: 14px;"><em><strong>Luca Brusa &middot; Catherine Matias</strong></em></span> </p>
<br>

We explore the use of hypergraphs to detect line clusters of points in $\mathbb{R}^2$. Similarly to the construction of pairwise similarity measures, we  here resort on third order affinity measures to detect alignment of points since pairwise measures would be useless to detect alignment.  Thus, for any triplet of points $\{i, j, k\}$, we use the mean distance to the best fitting line as a dissimilarity measure $d(i,j,k)$ and transform this through a Gaussian kernel to a similarity measure.


## Dataset generation 
There are 2 different experiments, generated through the Python Notebooks `1_synthetic_data_LINE_CLUSTERING_3groups.ipynb` (with 2 intersecting lines) and `1_synthetic_data_LINE_CLUSTERING_4groups.ipynb` (with 3 intersecting lines).

In each setting, we randomly generate the same number of points per line in the range $[-0.5, 0.5]^2$ and perturbed with centered Gaussian noise with standard deviation $0.01$. We then add noisy points, generated from uniform distribution on 
$[-0.5, 0.5]^2$.  
The particular settings of each  experiment is described in the following table

| Nb points per line | Noisy points | Total nb points |mean nb of hyperedges|
| --------|------|-----|------|
|2 lines | 30 | 40 | 100|  1070.84    |
| 3 lines | 30| 60 | 150|  587.7  |

For both settings, we randomly generated 100 3-uniform hypergraphs as follows. We repeatedly randomly drawn 3 points $\{i, j, k\}$, computed the mean distance $d(i,j,k)$  to the best fitting line 
and then their similarity using a Gaussian kernel $\exp(-d(i,j,k)^2/\sigma^2)$ with $\sigma^2=0.04$. We then construct a hyperedge $\{i,j,k\}$ whenever the similarity is larger than a  threshold $\epsilon= 0.999$. Note that in this way, we both construct signal hyperedges where all points come from same line cluster as well as noise hyperedges, where the points are sufficiently aligned without being issued from the same line. For each hypergraph, the signal:noise ratio of hyperedges, namely the ratio between signal and noise hyperedges is set to 2. We simulated sparse hypergraphs and the average number of hyperedges is shown in Table~\ref{tab:lines}.  Note that whenever there were isolated nodes in the hypergraph, we discarded them from the clustering analysis.  

## Clustering the nodes
We performed clustering of the nodes through 4 different methods: 

   - 2 modularity-based methods from Chodrow et al (2021) [1]: implemented in the Julia file `2_Chodrow_Line_Clustering.jl` that builds on the codes provided [on this GitHub page](https://github.com/nveldt/HyperModularity.jl)
   - the modularity-based method by Kaminski et al (2019) [2]: implemented in the Julia file `3_Kaminski_Line_Clustering.jl`  that builds on the codes provided [on this GitHub page](https://gist.github.com/pszufe/02666497d2c138d1b2de5b7f67784d2b)
   - the probabilistic model-based method HyperSBM from Matias and Brusa [3]: implemented in the R file `4_Line_Clustering_HSBM.R` that uses our R package available [on this GitHub page](https://github.com/LB1304/HyperSBM). This last file also contains the code to produce ARI plots and compare results. 
   

## References
  - [1] Chodrow et al (2021). Generative hypergraph clustering: From blockmodels to modularity. Science Advances, [DOI](https://doi.org/10.1126/sciadv.abh1303)
  - [2] Kaminski et al (2019). Clustering via hypergraph modularity. Plos One, [DOI](https://doi.org/10.1371/journal.pone.0224307)
  - [3] Brusa and Matias (2023). Model-based clustering in simple hypergraphs through a stochastic blockmodel. [ArXiV preprint](https://arxiv.org/abs/2210.05983)






