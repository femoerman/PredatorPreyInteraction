# t-SNE for predator traits

#0) Clear memory
rm(list=ls())

#1) Read in libraries
{
  library(ellipse)
  library(ggbiplot)
  library(reshape2)
  library(Rmisc)
  library(Rtsne)
  library(tidyverse)
  library(viridis)
}

#2) Read in and look at data
{
  lifehist = read.table("lifehist_summary.txt", header = T)
  dim(lifehist)
  colnames(lifehist) = c( "pred", "bact", "spec", "replicate", "trait", "value")
  traits = read.table("trait_summary.txt", header = T)
  traits
}

#3) Combine life history and trait data, and prepare data frame for analysis
{
  #3.1) Combine data
  alldata = rbind(lifehist, traits)
  head(alldata)
  dim(alldata)
  alldata = summarySE(na.omit(alldata), measurevar = "value", groupvars = c("pred", "bact", "spec", "trait", "replicate"))
  head(alldata)
  alldata = dcast(alldata, pred + bact + spec + replicate ~ trait, value.var = "value")
  head(alldata)
  
  #3.2) Replace NA values with column mean
  x = colnames(alldata[,-c(1:4)])
  x
  for(i in 1:length(x)){
    alldata[,x[i]][is.na(alldata[,x[i]])] = mean(alldata[,x[i]], na.rm = T)
  }
  alldata
  
  #3.3) Extract only data
  meta = alldata[,c(1:4)]
  data = alldata[,-c(1:4)]
  head(meta)
  head(data)
}

#4) Perform Rtsne function (may take some minutes to complete...)
{
  set.seed(9)
  tsne_model_1 = Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity=10, theta=0.5, dims=2)
  
  # Getting the two dimension matrix
  
  d_tsne_1 = as.data.frame(tsne_model_1$Y)
  
  d_tsne_1 = cbind(d_tsne_1, meta)
  head(d_tsne_1)
}

# #5) Plotting the results by adding treatment information
# {
#   mypal = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
#             "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
#             "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
#             "#8A7C64", "#599861")
#   
#   d_tsne_1$comb = paste(d_tsne_1$bact, d_tsne_1$pred, sep = "_")
#   
#   ggplot(d_tsne_1, aes(x=V1, y=V2, colour = pred)) +
#     geom_point(size=1) +
#     guides(colour=guide_legend(override.aes=list(size=6))) +
#     xlab("") + ylab("") +
#     ggtitle("t-SNE") +
#     scale_color_manual(values = mypal) + 
#     facet_grid(~spec) +
#     theme_bw() +
#     labs(colour = "Bacterial state") +
#     xlab("t-SNE dimension 1") +
#     ylab("t-SNE dimension 2") +
#     theme(legend.direction = "horizontal", 
#           legend.position = "bottom",
#           legend.box = "horizontal",
#           strip.background = element_rect(colour = "white", fill = "white"),
#           strip.text = element_text(face = "italic", colour = "black"),
#           panel.border = element_rect(colour = "grey"),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank())
#   
#   #ggsave("../../../manuscript/figures/t-SNE_pheno3_phage_allres.pdf")
# }

#6) Perform tsne separately for species
{
  
  #6.1) Bd
  {
    meta = alldata[alldata$spec == "Bd",c(1:4)]
    data = alldata[alldata$spec == "Bd",-c(1:4)]
    head(meta)
    head(data)
    
    # Rtsne function may take some minutes to complete...
    
    set.seed(9)
    tsne_model_1 = Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity=3, theta=0.5, dims=2)
    
    # Getting the two dimension matrix
    
    d_tsne_Bd = as.data.frame(tsne_model_1$Y)
    
    d_tsne_Bd = cbind(d_tsne_Bd, meta)
    head(d_tsne_Bd)
    
    # adding treatment information
    
    d_tsne_Bd$comb = paste(d_tsne_Bd$bact, d_tsne_Bd$pred, sep = "_")
    
  }
  
  #6.2) Ct
  {
    meta = alldata[alldata$spec == "Ct",c(1:4)]
    data = alldata[alldata$spec == "Ct",-c(1:4)]
    head(meta)
    head(data)
    
    # Rtsne function may take some minutes to complete...
    
    set.seed(9)
    tsne_model_1 = Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity=3, theta=0.5, dims=2)
    
    # Getting the two dimension matrix
    
    d_tsne_Ct = as.data.frame(tsne_model_1$Y)
    
    d_tsne_Ct = cbind(d_tsne_Ct, meta)
    head(d_tsne_Ct)
    
    # adding treatment information
    
    d_tsne_Ct$comb = paste(d_tsne_Ct$bact, d_tsne_Ct$pred, sep = "_")
    
  }
  
  #6.3) Ec
  {
    meta = alldata[alldata$spec == "Ec",c(1:4)]
    data = alldata[alldata$spec == "Ec",-c(1:4)]
    head(meta)
    head(data)
    
    # Rtsne function may take some minutes to complete...
    
    set.seed(9)
    tsne_model_1 = Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity=3, theta=0.5, dims=2)
    
    # Getting the two dimension matrix
    
    d_tsne_Ec = as.data.frame(tsne_model_1$Y)
    
    d_tsne_Ec = cbind(d_tsne_Ec, meta)
    head(d_tsne_Ec)
    
    # adding treatment information
    
    d_tsne_Ec$comb = paste(d_tsne_Ec$bact, d_tsne_Ec$pred, sep = "_")
    
  }
  
  #6.4) Jl
  {
    meta = alldata[alldata$spec == "Jl",c(1:4)]
    data = alldata[alldata$spec == "Jl",-c(1:4)]
    head(meta)
    head(data)
    
    # Rtsne function may take some minutes to complete...
    
    set.seed(9)
    tsne_model_1 = Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity = 3, theta=0.5, dims=2)
    
    # Getting the two dimension matrix
    
    d_tsne_Jl = as.data.frame(tsne_model_1$Y)
    
    d_tsne_Jl = cbind(d_tsne_Jl, meta)
    head(d_tsne_Jl)
    
    # adding treatment information
    
    d_tsne_Jl$comb = paste(d_tsne_Jl$bact, d_tsne_Jl$pred, sep = "_")
    
  }
  
  #6.5) Pf
  {
    meta = alldata[alldata$spec == "Pf",c(1:4)]
    data = alldata[alldata$spec == "Pf",-c(1:4)]
    head(meta)
    head(data)
    
    # Rtsne function may take some minutes to complete...
    
    set.seed(9)
    tsne_model_1 = Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity=3, theta=0.5, dims=2)
    
    # Getting the two dimension matrix
    
    d_tsne_Pf = as.data.frame(tsne_model_1$Y)
    
    d_tsne_Pf = cbind(d_tsne_Pf, meta)
    head(d_tsne_Pf)
    
    # adding treatment information
    
    mypal = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
              "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
              "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
              "#8A7C64", "#599861")
    
    d_tsne_Pf$comb = paste(d_tsne_Pf$bact, d_tsne_Pf$pred, sep = "_")
    
  }
  
  #6.6) Sc
  {
    meta = alldata[alldata$spec == "Sc",c(1:4)]
    data = alldata[alldata$spec == "Sc",-c(1:4)]
    head(meta)
    head(data)
    
    # Rtsne function may take some minutes to complete...
    
    set.seed(9)
    tsne_model_1 = Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity=3, theta=0.5, dims=2)
    
    # Getting the two dimension matrix
    
    d_tsne_Sc = as.data.frame(tsne_model_1$Y)
    
    d_tsne_Sc = cbind(d_tsne_Sc, meta)
    head(d_tsne_Sc)
    
    # adding treatment information
    
    d_tsne_Sc$comb = paste(d_tsne_Sc$bact, d_tsne_Sc$pred, sep = "_")
  }
  
  #6.7 Sm
  {
    meta = alldata[alldata$spec == "Sm",c(1:4)]
    data = alldata[alldata$spec == "Sm",-c(1:4)]
    head(meta)
    head(data)
    
    # Rtsne function may take some minutes to complete...
    
    set.seed(9)
    tsne_model_1 = Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity=3, theta=0.5, dims=2)
    
    # Getting the two dimension matrix
    
    d_tsne_Sm = as.data.frame(tsne_model_1$Y)
    
    d_tsne_Sm = cbind(d_tsne_Sm, meta)
    head(d_tsne_Sm)
    
    # adding treatment information
    
    d_tsne_Sm$comb = paste(d_tsne_Sm$bact, d_tsne_Sm$pred, sep = "_")
    
  }
  
  #6.8) Combine the results of all analyses
  {
    tsne_all = rbind(d_tsne_Bd,
                     d_tsne_Ct,
                     d_tsne_Ec,
                     d_tsne_Jl,
                     d_tsne_Pf,
                     d_tsne_Sc,
                     d_tsne_Sm)
    
    tsne_all$comb = factor(tsne_all$comb, levels = c("A_A", "E_A", "A_E", "E_E"),
                           labels = c("Ancestral predator, ancestral prey",
                                      "Ancestral predator, evolved prey",
                                      "Evolved predator, ancestral prey",
                                      "Evolved predator, Evolved prey"))
    
    tsne_all$spec = factor(tsne_all$spec, levels = c("Bd", "Ct", "Ec", "Jl", "Pf", "Sc", "Sm"),
                           labels = c("Brevundimonas", "Comamonas", "Escherichia", "Janthinobacterium", "Pseudomonas", "Sphingomonas", "Serratia"))
  }
  
  #6.9) Plot the combined results
  {
  ggplot(tsne_all, aes(x=V1, y=V2)) +
    #geom_point(shape=21, alpha = 0.8, size=5, colour = "black") +
    geom_point(aes(fill=pred,shape=bact), alpha = 0.7, size=5, colour = "black")+
    scale_shape_manual(values=c(21, 22), labels = c("Ancestral", "Evolved"))+
    xlab("t-SNE dimension 1") +
    ylab("t-SNE dimension 2") +
    #scale_fill_viridis(discrete = T, direction = -1, option = "plasma", end = 0.8) +
    #scale_colour_viridis(discrete = T, direction = -1, option = "plasma", end = 0.8) +
    scale_fill_manual(values = c("#016392", "#c72e29"), labels = c("Ancestral", "Evolved")) +
    #scale_colour_manual(values = c("#016392", "#c72e29")) +
    theme_bw() +
    facet_wrap(~spec, scales = "free", ncol = 4) +
    labs(fill = "Predator evolutionary\nhistory",shape = "Prey evolutionary\nhistory") +
    guides(fill= guide_legend(override.aes=list(shape=c(21)), order = 1)) +
    theme(axis.text.x=element_text(size = 14),
          axis.text.y=element_text(size = 14),
          axis.title.x=element_text(size = 20),
          axis.title.y=element_text(size = 20),
          strip.background = element_rect(colour = "white", fill = "white"),
          strip.text = element_text(face = "italic", colour = "black", size = 20),
          panel.border = element_rect(colour = "grey", size = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(0.93, 0.2),
          legend.justification = c(0.75, 0.4),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 20))
  
  ggsave("Figure_1.pdf", width = 15, height = 8)
  }
}