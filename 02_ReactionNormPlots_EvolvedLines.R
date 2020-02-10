#1) Clear memory
rm(list=ls())

#2) Load libraries
{
  library(tidyverse)
  library(ggthemes)
  library(lme4)
  library(nlme)
  library(MuMIn)
  library(reshape2)
  library(Rmisc)
}

#3) Set working directory and load dataset
{
  load("Posteriors_BacteriaData.RData")
  dd <- output$sumoutput
  
  #Combine posterior output with info on strains
  load("Population_Data_t0.RData")
  pop_output <- filter(pop_output, medium== "5% KB")
  
  dd <- cbind(dd, pop_output[, 2:6])
  
  #Create variables pred (e/a), bact(e/a) and spec (bact. sp.)
  dd$pred <- c(rep("A", 48), rep("E", 48))
  dd$bact <- c(rep("A", 24), rep("E", 24), rep("A", 24), rep("E", 24))
  dd$spec <- rep(c(rep("Ec", 3), rep("Jl", 3), rep("Sc", 3), rep("Bd", 3), rep("Pf", 3), rep("Ct", 6), rep("Sm", 3)), 4)
  dd$spec2 <- rep(c(rep("Ec", 3), rep("Jl", 3), rep("Sc", 3), rep("Bd", 3), rep("Pf", 3), rep("F", 3), rep("Ct", 3), rep("Sm", 3)), 4)
  

}


#4) Prepare data for final plotting

dd2 = melt(dd[,c("logr0.mean", "logK.mean", "logalpha.mean", "pred", "bact", "spec",  "spec2", "replicate")], id.vars = c("pred", "bact", "spec", "spec2", "replicate"))
colnames(dd2) = c("pred", "bact", "spec", "spec2", "replicate", "trait", "value")


#5) Do the linear models for r0, K and alpha
{
  #5.1) r0
  r0model <- lm(data=dd, formula = logr0.mean ~ bact*pred*spec)
  bestr0 <- step(r0model, direction = "backward")
  anova(bestr0)
  summary(bestr0)

  #5.2) alpha
  alphamodel <- lm(data=dd, formula = logalpha.mean ~ bact*pred*spec)
  bestalpha <- step(alphamodel, direction = "backward")
  anova(bestalpha)
  summary(bestalpha)

  #5.3) K
  Kmodel <- lm(data=dd, formula = logK.mean ~ bact*pred*spec)
  bestK <- step(Kmodel, direction = "backward")
  anova(bestK)
  summary(bestK)
}

#6) Extract predictions and confidence intervals from the models for plotting
{
  dd3 = summarySE(na.omit(dd2), measurevar = "value", groupvars = c("pred", "bact", "spec", "trait"))
  head(dd3)
  dd3_r = dd3[dd3$trait == "logr0.mean",]
  dd3_rp = predict(bestr0, newdata = dd3_r, interval = 'confidence')
  dd3_r = cbind(dd3_r, dd3_rp)
  dd3_r
  dd3_K = dd3[dd3$trait == "logK.mean",]
  dd3_Kp = predict(bestK, newdata = dd3_K, interval = 'confidence')
  dd3_K = cbind(dd3_K, dd3_Kp)
  dd3_K
  dd3_alpha = dd3[dd3$trait == "logalpha.mean",]
  dd3_alphap = predict(bestalpha, newdata = dd3_alpha, interval = 'confidence')
  dd3_alpha = cbind(dd3_alpha, dd3_alphap)
  dd3_alpha
  
  dd3 = rbind(dd3_r, dd3_K, dd3_alpha)
  head(dd3)
}

#7) Edit/prepare tables for nice plotting
{

  #7.1) dd2
  
  dd2 = na.omit(dd2[,c("pred", "bact", "spec", "trait", "value")])
  head(dd2)
  
  dd2$trait = gsub("logr.*", "italic(r[0])", dd2$trait)
  dd2$trait = gsub("logK.*", "italic(K)", dd2$trait)
  dd2$trait = gsub("logalpha.*", "alpha", dd2$trait)
  dd2$bact = as.numeric(as.character(ifelse(dd2$bact == "A", 1, 2)))
  
  dd2$trait = factor(dd2$trait, levels = c("italic(r[0])", "italic(K)", "alpha"),
                     labels = c(expression(paste("Intrinsic rate of increase, log ", italic(r[0]), " ", bgroup("(", frac(1,h), ")"))), expression(paste("Population equilibrium density, log ", italic(K), " ", bgroup("(", frac(individuals,mL),")"))), expression(paste("Competitive ability, log ", alpha, " ", bgroup("(", frac(mL,h%.%individuals),")")))))
  
  
  dd2$pred = ifelse(dd2$pred == "A", "Ancestral", "Evolved")
  
  dd2$spec = gsub("Bd", "italic(Brevundimonas)", dd2$spec)
  dd2$spec = gsub("Ct", "italic(Comamonas)", dd2$spec)
  dd2$spec = gsub("Ec", "italic(Escherichia)", dd2$spec)
  dd2$spec = gsub("Jl", "italic(Janthinobacterium)", dd2$spec)
  dd2$spec = gsub("Pf", "italic(Pseudomonas)", dd2$spec)
  dd2$spec = gsub("Sc", "italic(Sphingomonas)", dd2$spec)
  dd2$spec = gsub("Sm", "italic(Serratia)", dd2$spec)
  
  #7.2) dd3
  
  dd3$trait = gsub("logr.*", "italic(r[0])", dd3$trait)
  dd3$trait = gsub("logK.*", "italic(K)", dd3$trait)
  dd3$trait = gsub("logalpha.*", "alpha", dd3$trait)
  dd3$bact = as.numeric(as.character(ifelse(dd3$bact == "A", 1, 2)))
  
  dd3$trait = factor(dd3$trait, levels = c("italic(r[0])", "italic(K)", "alpha"),
                     labels = c(expression(paste("Intrinsic rate of increase, log ", italic(r[0]), " ", bgroup("(", frac(1,h), ")"))), expression(paste("Population equilibrium density, log ", italic(K), " ", bgroup("(", frac(individuals,mL),")"))), expression(paste("Competitive ability, log ", alpha, " ", bgroup("(", frac(mL,h%.%individuals),")")))))
  
  dd3$pred = ifelse(dd3$pred == "A", "Ancestral", "Evolved")
  
  dd3$spec = gsub("Bd", "italic(Brevundimonas)", dd3$spec)
  dd3$spec = gsub("Ct", "italic(Comamonas)", dd3$spec)
  dd3$spec = gsub("Ec", "italic(Escherichia)", dd3$spec)
  dd3$spec = gsub("Jl", "italic(Janthinobacterium)", dd3$spec)
  dd3$spec = gsub("Pf", "italic(Pseudomonas)", dd3$spec)
  dd3$spec = gsub("Sc", "italic(Sphingomonas)", dd3$spec)
  dd3$spec = gsub("Sm", "italic(Serratia)", dd3$spec)
}

# 8) Plot Model predictions along with raw data
{
  dd3
  p1 = ggplot(dd3, aes(bact, fit, colour = pred, fill = pred)) +
    geom_line(size = 2) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, colour = NULL), alpha = 0.5) +
    geom_point(dd2, mapping = aes(bact, value), position = position_dodge(width = 0.15), size = 3, alpha = 0.5) +
    facet_grid(trait~spec, scale = "free", labeller= label_parsed, switch = "y") +
    theme_classic()+ scale_colour_manual(values = c("#016392", "#c72e29")) + scale_fill_manual(values = c("#016392", "#c72e29")) +
    scale_x_continuous(breaks = c(1, 2), labels = c("Ancestral", "Evolved"), expand = c(0.1, 0.1)) +
    labs(colour = "Evolutionary\nhistory of\npredator",
         fill = "Evolutionary\nhistory of\npredator") +
    xlab("Evolutionary history of prey") +
    ylab("Life history traits") +
    theme(panel.spacing = unit(2, "lines"),
          axis.title=element_text(size=14),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 20),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), size = 20),
          strip.background = element_rect(colour = "white", fill = "white"),
          strip.text.x = element_text(face = "italic", colour = "black", size = 15),
          strip.text.y = element_text(colour = "black", size = 15),
          strip.placement = "outside",
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 20))#,
          #panel.grid.minor.y = element_line(color = "grey80", linetype = 2))
  p1
}
ggsave("Figure_2.pdf", width = 20, height = 16)
