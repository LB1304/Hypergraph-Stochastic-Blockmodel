###########
# HSBM on LineClustering simulations
###########

#require(devtools)
#devtools::install_github("LB1304/HyperSBM")
require(HyperSBM)
REP= 99 


# 2linecluster
setwd("data_2linecluster/")
for (rep in 0:REP){
  # import true clusters 
  filename = paste0("rep_", rep,"/true_clusters.txt")
  classes <- paste(readLines(filename), collapse =",")
  Ztrue_2_rep <- as.numeric(unlist(strsplit(classes, ",")))
  
  ## Import dataset 
  HG <- HyperSBM::import_Hypergraph(file_name = paste0("rep_",rep,"/hyperedges.txt"))
  ## Estimate with full model (model=0) - initialization with ssc (start=2)
  res <- vector(mode = "list", length = 6)
  ICL <- c()
  for (q in 1:6){
    print(paste0("repetition ", rep, " with q=",q))
    res[[q]]  <- HyperSBM::HSBM(Hypergraph = HG, Q = q, start = 2, model = 0, tol = 1e-5, maxit_VEM = 100, maxit_FP = 25, n_threads = 25, print = TRUE)
    ICL <- c(ICL,res[[q]]$ICL)
  }
  
  ## Select best model and compute ARI
  Qhat <-which.max(ICL)
  ariHSBM <- ari_index(res[[Qhat]]$Z, Ztrue_2_rep)
  
  ## save results 
  save(res, ICL, Qhat, ariHSBM, file = paste0("rep_", rep, "_HSBMres_Qsel.RData"))
  print(paste0("repetition ", rep, " finished"))
}


# 3linecluster
setwd("data_3linecluster/")
for (rep in 0:REP){

  ## import true clusters 
  filename = paste0("rep_", rep,"/true_clusters.txt")
  classes <- paste(readLines(filename), collapse =",")
  Ztrue_3_rep <- as.numeric(unlist(strsplit(classes, ",")))
  
  ## Import dataset 
  HG <- HyperSBM::import_Hypergraph(file_name = paste0("rep_",rep,"/hyperedges.txt") )
  
  ## Estimate with full model (model=0) - initialization with ssc (start=2)
  res <- vector(mode = "list", length = 6)
  ICL <- c()
  for (q in 1:6){
    print(paste0("repetition ", rep, " with q=",q))
    res[[q]]  <- HyperSBM::HSBM(Hypergraph = HG, Q = q, start = 2, model = 0, tol = 1e-5, maxit_VEM = 100, maxit_FP = 25, n_threads = 25, print = TRUE)
    ICL <- c(ICL,res[[q]]$ICL)
  }
  
  ## Select best model and compute ARI
  Qhat <-which.max(ICL)
  ariHSBM <- ari_index(res[[Qhat]]$Z, Ztrue_3_rep)
  
  ## save results 
  save(res, ICL, Qhat, ariHSBM, file = paste0("rep_", rep, "_HSBMres_Qsel.RData"))
  print(paste0("repetition ", rep, " finished"))
}


############
# Plots
############
library(ggplot2)


#####
## 2lineclustering 
#####

# Chodrow's results 
Chod_AON <- as.vector(read.table("Chodrow_AON_ARI_2lines.txt"))
Chod_Symm <- as.vector(read.table("Chodrow_Symm_ARI_2lines.txt"))

# Kaminski's results 
Kam <- as.vector(read.table("Kaminski_ARI_2lines.txt"))

# HSBM results  
REP=99
setwd("data2_HSBMres_Qsel")
ariHSBM_global<-c()
Qhat_global <- c()
ICL_global <- c()
for (rep in 0:REP){
  load(paste0("rep_", rep, "_HSBMres_Qsel.RData"))
  ariHSBM_global <- c(ariHSBM_global,ariHSBM)
  ICL_global <- rbind(ICL_global,ICL) 
  Qhat_global <- c(Qhat_global,Qhat)
}

### ARI Plot
ARI <- data.frame(ari=as.numeric(c(Chod_AON$V1,Chod_Symm$V1,Kam$V1,ariHSBM_global)),method=c(rep("Chodrow_AON",100),rep("Chodrow_Symm",100),rep("Kaminski",100),rep("HSBM",100)))
png(file = "ARI_2lines.png")
ggplot(data=ARI)+aes(x=method,y=ari) +geom_boxplot()
dev.off()

png(file = "ARI_2lines_color.png")
ggplot(data=ARI)+aes(x=method,y=ari,col=method) +geom_boxplot()
dev.off()


### ICL plot of all curves - Not Nice
ICL_df <-  data.frame(Q = rep(1:6, times = 100), ICL = as.vector(t(ICL_global)), replicate=rep(1:100, each = 6))
ggplot(data=ICL_df, aes(x=Q,y=ICL,group=replicate)) +geom_line() +theme_bw()


#### Look at the number of groups estimated by  methods 
setwd("..")
Chod_AON_Q <- as.vector(read.table("Chodrow_AON_Qhat_2lines.txt"))
Chod_Symm_Q <- as.vector(read.table("Chodrow_Symm_Qhat_2lines.txt"))
Kam_Q <- as.vector(read.table("Kaminski_Qhat_2lines.txt"))

Qhat <- data.frame(Q=as.factor(c(Chod_AON_Q$V1,Chod_Symm_Q$V1,Kam_Q$V1,Qhat_global)),method=c(rep("Chodrow_AON",100),rep("Chodrow_Symm",100),rep("Kaminski",100),rep("HSBM",100)))
png(file = "Qhat_2lines.png")
ggplot(data=Qhat ) + 
  geom_bar(aes(x = Q, group=method, fill=method), position = "dodge" )+
  scale_fill_grey()+  
  theme_bw()+
  theme(legend.position = c(1, 1),
           legend.justification = c(1, 1))
dev.off()

png(file = "Qhat_2lines_color.png")
ggplot(data=Qhat ) + 
  geom_bar(aes(x = Q, group=method, fill=method), position = "dodge" )+
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1))
dev.off()



#####
## 3lineclustering 
#####

# Chodrow's results 
Chod_AON <- as.vector(read.table("Chodrow_AON_ARI_3lines.txt"))
Chod_Symm <- as.vector(read.table("Chodrow_Symm_ARI_3lines.txt"))

# Kaminski's results 
Kam <- as.vector(read.table("Kaminski_ARI_3lines.txt"))

# HSBM results  
REP=99
setwd("data3_HSBMres_Qsel")
ariHSBM_global<-c()
Qhat_global <- c()
ICL_global <- c()
for (rep in 0:REP){
  load(paste0("rep_", rep, "_HSBMres_Qsel.RData"))
  ariHSBM_global <- c(ariHSBM_global,ariHSBM)
  ICL_global <- rbind(ICL_global,ICL) 
  Qhat_global <- c(Qhat_global,Qhat)
}

### ARI Plot
ARI <- data.frame(ari=as.numeric(c(Chod_AON$V1,Chod_Symm$V1,Kam$V1,ariHSBM_global)),method=c(rep("Chodrow_AON",100),rep("Chodrow_Symm",100),rep("Kaminski",100),rep("HSBM",100)))
png(file = "ARI_3lines.png")
ggplot(data=ARI)+aes(x=method,y=ari) +geom_boxplot()
dev.off()

png(file = "ARI_3lines_color.png")
ggplot(data=ARI)+aes(x=method,y=ari,col=method) +geom_boxplot()
dev.off()


### ICL plot of all curves - Not Nice
ICL_df <-  data.frame(Q = rep(1:6, times = 100), ICL = as.vector(t(ICL_global)), replicate=rep(1:100, each = 6))
ggplot(data=ICL_df, aes(x=Q,y=ICL,group=replicate)) +geom_line() +theme_bw()


#### Look at the number of groups estimated by  methods 
setwd("..")
Chod_AON_Q <- as.vector(read.table("Chodrow_AON_Qhat_3lines.txt"))
Chod_Symm_Q <- as.vector(read.table("Chodrow_Symm_Qhat_3lines.txt"))
Kam_Q <- as.vector(read.table("Kaminski_Qhat_3lines.txt"))

Qhat <- data.frame(Q=as.factor(c(Chod_AON_Q$V1,Chod_Symm_Q$V1,Kam_Q$V1,Qhat_global)),method=c(rep("Chodrow_AON",100),rep("Chodrow_Symm",100),rep("Kaminski",100),rep("HSBM",100)))


png(file = "Qhat_3lines.png")
ggplot(data=Qhat ) + 
  geom_bar(aes(x = Q, group=method, fill=method), position = "dodge" )+
  scale_fill_grey()+  
  theme_bw()+
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1))
dev.off()

png(file = "Qhat_3lines_color.png")
ggplot(data=Qhat ) + 
  geom_bar(aes(x = Q, group=method, fill=method), position = "dodge" )+
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1))
dev.off()

