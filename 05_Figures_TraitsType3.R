#1) Clear memory
rm(list=ls())

#2) Load libraries
{
  library(tidyverse)
  library(MuMIn)
  library(car)
  library(Rmisc)
}

#3) Set working directory and load dataset
{
  
  #3.1) Load the datasets
  {
    load("rawbactdata.RData")                              #Raw protist dataset
    dd<- databact
    load("Posteriors_BacteriaData.RData")            #Posterior data from model fits
    sumdata <- output$sumoutput
  }
  
  #3.2) Create variables pred (e/a), bact(e/a) and spec (bact. sp.)
  {
    dd$pred <- c(rep("A", 48), rep("E", 48))
    dd$bact <- c(rep("A", 24), rep("E", 24), rep("A", 24), rep("E", 24))
    dd$spec <- rep(c(rep("Ec", 3), rep("Jl", 3), rep("Sc", 3), rep("Bd", 3), rep("Pf", 3), rep("Ct", 6), rep("Sm", 3)), 4)
    dd$spec2 <- rep(c(rep("Ec", 3), rep("Jl", 3), rep("Sc", 3), rep("Bd", 3), rep("Pf", 3), rep("F", 3), rep("Ct", 3), rep("Sm", 3)), 4)
    dd$indiv_per_ml <- dd$indiv_per_volume*1000
    sumdata$pred <- c(rep("A", 48), rep("E", 48))
    sumdata$bact <- c(rep("A", 24), rep("E", 24), rep("A", 24), rep("E", 24))
    sumdata$spec <- rep(c(rep("Ec", 3), rep("Jl", 3), rep("Sc", 3), rep("Bd", 3), rep("Pf", 3), rep("Ct", 6), rep("Sm", 3)), 4)
    sumdata$spec2 <- rep(c(rep("Ec", 3), rep("Jl", 3), rep("Sc", 3), rep("Bd", 3), rep("Pf", 3), rep("F", 3), rep("Ct", 3), rep("Sm", 3)), 4)
  }
  
  #3.3) Import bacterial cell count data (Flow cytometry) and concatenate to the rest of the data
  {
    #Import bacterial data
    bactdata <- read.csv("bacteria_data_SH_20170420.csv", sep=",")
    
    #Remove the data for the axenic cultures and those of the last timepoint, where we don't have protist data
    bact2 <- subset(bactdata, Bactsp!="-")
    bact2 <- filter(bact2, measument != "t9")
    
    #paste bacterial density data to rest of the data
    dd$bactcounts <- bact2$counts_ml
  }
  
  #3.4) Create data sets for t=1, mid-log phase and 99.9%K
  {
    #Make dataset for t1
    dd.t0 <- filter(dd, measument=="t1")
    
    #Calculate mid log phase and K for all populations
    {
      sumdata <- drop_na(sumdata)
      dd.K <- data.frame()
      dd.inf <- data.frame()
      
      
      row.names(sumdata) <- sumdata$ident
      for (i in sumdata$ident){
        tempdata <- filter(dd, ident==i)
        timeK <- min(filter(tempdata, indiv_per_ml>=0.99*exp(sumdata[i, ]$logK.mean))$time)
        timeK <- ifelse(timeK==Inf, tempdata[which.max(tempdata$indiv_per_ml), ]$time, timeK)
        if(timeK> 0){
          timeInflect <- tempdata[which.min((tempdata[which(tempdata$time<timeK), ]$indiv_per_ml - exp(sumdata[i, ]$logK.mean)/2)^2), ]$time
          dd.inf <- rbind(dd.inf, filter(tempdata, time==timeInflect))
        }
        dd.K <- rbind(dd.K, filter(tempdata, time==timeK))
      }
    }
  }
  
  #3.5) remove datapoints where no good datapoint was available (i.e. where population crashed and never reached K or mid log phase)
  dd.K <- na.omit(dd.K)
  dd.inf <- na.omit(dd.inf)
}

#4) Fit models over all datapoints with Bacterial cell counts and protist cell counts as covariates
{
  #4.1)major cell axis size
  {
     #Model selection
    tall.major.mean.full <- lm(data=dd, formula = major_mean ~ bact*pred*spec*log(indiv_per_ml)*log(bactcounts))
    tall.major.mean.null <- lm(data=dd, formula = major_mean ~ 1)
    besttall.major.mean1 <- step(object = tall.major.mean.null, scope = list(upper=tall.major.mean.full, lower = tall.major.mean.null), direction = "both")
    besttall.major.mean2 <- step(object = tall.major.mean.full, scope = list(upper=tall.major.mean.full, lower = tall.major.mean.null), direction = "both")
    AICc(besttall.major.mean1, besttall.major.mean2)
    summary(besttall.major.mean1)
    Anova(besttall.major.mean1, type="III", contrasts = c("contr.sum", "contr.poly"))
  }
  
  #4.2) gross speed of cells
  {
    #Model selection
    tall.gross.speed.full <- lm(data=dd, formula = gross_speed_mean ~ bact*pred*spec*log(indiv_per_ml)*log(bactcounts))
    tall.gross.speed.null <- lm(data=dd, formula = gross_speed_mean ~ 1)
    besttall.gross.speed1 <- step(object = tall.gross.speed.null, scope = list(upper=tall.gross.speed.full, lower = tall.gross.speed.null), direction = "both")
    besttall.gross.speed2 <- step(object = tall.gross.speed.full, scope = list(upper=tall.gross.speed.full, lower = tall.gross.speed.null), direction = "both")
    AICc(besttall.gross.speed1, besttall.gross.speed2)
    Anova(besttall.gross.speed1, type="III", contrasts = c("contr.sum", "contr.poly"))
    summary(besttall.gross.speed1)
  }
  
  #4.3) Turning angle distribution of cells
  {
   #Model selection
    tall.sd.turning.full <- lm(data=dd, formula = log10(sd_turning_mean) ~ bact*pred*spec*log(indiv_per_ml)*log(bactcounts))
    tall.sd.turning.null <- lm(data=dd, formula = log10(sd_turning_mean) ~ 1)
    besttall.sd.turning1<- step(object = tall.sd.turning.null, scope = list(upper=tall.sd.turning.full, lower = tall.sd.turning.null), direction = "both")
    besttall.sd.turning2<- step(object = tall.sd.turning.full, scope = list(upper=tall.sd.turning.full, lower = tall.sd.turning.null), direction = "both")
    AICc(besttall.sd.turning1, besttall.sd.turning2)
    Anova(besttall.sd.turning1, type="III", contrasts = c("contr.sum", "contr.poly"))
    summary(besttall.sd.turning1)
  }
}


#5) Prepare prediction data for plotting
{
  #Add all data together
  dd.t0$level <- "Exponential phase"
  dd.inf$level <- "Inflection point"
  dd.K$level <- "Equilibrium population density"
  dd.all <- rbind(dd.t0, dd.inf, dd.K)
  
  #Create prediction data
  predictions <-expand.grid(pred = unique(dd.all$pred), bact = unique(dd.all$bact), spec = unique(dd.all$spec), 
                        indiv_per_ml = c(mean(dd.t0$indiv_per_ml), mean(dd.inf$indiv_per_ml), mean(dd.K$indiv_per_ml)),
                        bactcounts = c(quantile(dd.all$bactcounts, 0.05), quantile(dd.all$bactcounts, 0.5), quantile(dd.all$bactcounts, 0.95)))
  
  #_Make predictions
  prd_size <- predict(besttall.major.mean1, newdata = predictions, se.fit = T)
  prd_speed <- predict(besttall.gross.speed1, newdata = predictions, se.fit = T)
  prd_turning <- predict(besttall.sd.turning1, newdata = predictions, se.fit = T)
  
  #Concatenate with prediction variables
  {
    #For size
    predictions$meansize = prd_size$fit
    predictions$uppersize = prd_size$fit + 1.96*prd_size$se.fit
    predictions$lowersize = prd_size$fit - 1.96*prd_size$se.fit
    
    #For speed
    predictions$meanspeed = prd_speed$fit
    predictions$upperspeed = prd_speed$fit + 1.96*prd_speed$se.fit
    predictions$lowerspeed = prd_speed$fit - 1.96*prd_speed$se.fit
    
    #For turning
    predictions$meanturning = prd_turning$fit
    predictions$upperturning = prd_turning$fit + 1.96*prd_turning$se.fit
    predictions$lowerturning = prd_turning$fit - 1.96*prd_turning$se.fit
  }
  
  #Add timepoint info to predictions
  predictions$level <- ifelse(predictions$indiv_per_ml == min(predictions$indiv_per_ml), "Exponential phase", 
                              ifelse(predictions$indiv_per_ml==max(predictions$indiv_per_ml), "Equilibrium population density", "Inflection point"))
  
}


#6) Create predictions for average prey and predator density
{
  prd.data <- expand.grid(bactcounts = quantile(dd$bactcounts, c(0.05,0.5, 0.95)), indiv_per_ml = quantile(dd$indiv_per_ml[dd$indiv_per_ml>0], c(0.05,0.5, 0.95)), spec = unique(dd$spec),
                          bact = unique(dd$bact), pred = unique(dd$pred))
  prd.data$size = predict(besttall.major.mean1, newdata = prd.data)
  prd.data$speed = predict(besttall.gross.speed1, newdata = prd.data)
  prd.data$turning = exp(predict(besttall.sd.turning1, newdata = prd.data))
  detach(package:Rmisc)
  detach(package:plyr)
  sum <- prd.data %>%
    group_by(spec, bact, bactcounts, indiv_per_ml) %>%
    summarise(size_ratio=size[2]/size[1], speed_ratio = speed[2]/speed[1], turning_ratio = turning[2]/turning[1])
  sum2 <- mutate(sum, "Log (Prey density)" = log(bactcounts), "Log (Predator density)" = log(indiv_per_ml))
  write.csv2(sum2, file = "PredictionsMorphologicalTraits.csv")

  sum3 <- gather(sum2, key = "Trait", value = "ratio", size_ratio, speed_ratio, turning_ratio)
  sum3$Trait <- ifelse(sum3$Trait=="size_ratio", "Cell size", ifelse(sum3$Trait=="speed_ratio", "Gross cell speed", "Turning angle"))
  sum3$'Predatordensity' <- ifelse(sum3$'Log (Predator density)' == min(sum3$'Log (Predator density)'), "5 %", ifelse(sum3$'Log (Predator density)' == max(sum3$'Log (Predator density)'), "95 %", "50 %"))
  sum3$'Preydensity' <- ifelse(sum3$'Log (Prey density)' == min(sum3$'Log (Prey density)'), "5 % prey density", ifelse(sum3$'Log (Prey density)' == max(sum3$'Log (Prey density)'), "95 % prey density", "50 % prey density"))
}

#7) Plot predictions (as ratios of the predicted trait values between ancestral and evolved predators)

ggplot(sum3, aes(x = Predatordensity, y = ratio, colour = bact))  + 
  geom_point(size=3, alpha = 0.5, position = position_dodge(width = 0.15)) + 
  facet_grid(Trait~Preydensity, switch = "y") +
  xlab("Predator density (quantile)") +
  ylab(expression("Trait ratio ("*frac("Evolved predator", "Ancestral predator")*")")) +
  theme_classic() + 
  scale_colour_manual(values = c("#016392", "#c72e29"), 
                      labels = c("Ancestral", "Evolved"), 
                      name = "Evolutionary history of prey") +
  theme(panel.spacing = unit(2, "lines"),
        axis.title=element_text(size=20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(colour = "black", size = 16),
        strip.placement = "outside",
        legend.position = "bottom", 
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 20))#,
        #panel.grid.minor.y = element_line(color = "grey80", linetype = 2))

ggsave("Figure_3.pdf", width = 15, height = 10)
