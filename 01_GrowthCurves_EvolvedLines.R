#1) Clear memory
{
  rm(list=ls())
}

#2) Load libraries
{
  library(tidyverse)
  # library(devtools)
  # install_github("femoerman/PBPGM")
  library(PBPGM)
}

#3) Set working directory and read in data
{
  load("2_data/GrowthCurves/Population_Data_t0.RData")
  dd <- pop_output
  load("2_data/GrowthCurves/Population_Data_t1.RData")
  dd <- rbind(dd, pop_output)
  load("2_data/GrowthCurves/Population_Data_t2.RData")
  dd <- rbind(dd, pop_output)
  load("2_data/GrowthCurves/Population_Data_t3.RData")
  dd <- rbind(dd, pop_output)
  load("2_data/GrowthCurves/Population_Data_t4.RData")
  dd <- rbind(dd, pop_output)
  load("2_data/GrowthCurves/Population_Data_t5.RData")
  dd <- rbind(dd, pop_output)
  load("2_data/GrowthCurves/Population_Data_t6.RData")
  dd <- rbind(dd, pop_output)
  load("2_data/GrowthCurves/Population_Data_t7.RData")
  dd <- rbind(dd, pop_output)
  load("2_data/GrowthCurves/Population_Data_t8.RData")
  dd <- rbind(dd, pop_output)
}

#4) Determine number of hours since start and indiv_per_ml
{
  dd$hours <- (as.numeric(as.Date(dd$date, "%Y-%m-%d"))-17266)*24 + dd$time-dd$time[1:125]
  dd$indiv_per_ml <- dd$indiv_per_volume*1000
}

#5) Plot the raw data
{
  ggplot(dd, aes(x=hours, y=indiv_per_ml, group=file)) + geom_point() + geom_line() + facet_wrap(~medium)
  
  #Add colour based on treatment
  ggplot(filter(dd, medium=="5% KB"), aes(x=hours, y=indiv_per_ml, group=file, colour=bacteria)) + geom_point() + geom_line() + facet_wrap(~predator)

  #Plot by bacteria and add colour based on Tetrahymena
  ggplot(filter(dd, medium=="5% KB"), aes(x=hours, y=indiv_per_ml, group=file, colour=predator)) + geom_point() + geom_line() + facet_wrap(~bacteria)
  
}

#6) Prepare data to fit growthcurves of populations with bacteria
{
  databact <- filter(dd, medium=="5% KB") %>% rename(clocktime=time)
  databact <- rename(databact, time=hours, popsize=indiv_per_ml, ident=file)
  #save bacterial data
  save(databact, file="2_data/rawbactdata.RData")
}

#7) Fit Beverton Holt model for all the population growth curves
{
  output <- BevertonHolt(databact, K.prior=log(1e4), Ksd.prior = 0.5, r0.prior = -2.3, r0sd.prior = 0.5, N0.prior = 2.3, 
                         N0sd.prior = 1, d.prior = -2.3, dsd.prior = 1, outputtype = "both", graphname = "2_data/ModelFitBacteria.pdf")
  
  
  save(output, file="2_data/Posteriors_BacteriaData.RData")
}

#8) Prepare data to fit growthcurves of populations without bacteria
{
  datacont <- filter(dd, medium!="5% KB") %>% rename(clocktime=time)
  datacont <- rename(datacont, time=hours, popsize=indiv_per_ml, ident=file)
}

#7) Fit Beverton Holt model for all the population growth curves
{
  output2 <- BevertonHolt(datacont, K.prior=log(7e4), Ksd.prior = 0.5, r0.prior = -2.3, r0sd.prior = 0.5, N0.prior = 2.3, 
                         N0sd.prior = 1, d.prior = -2.3, dsd.prior = 1, outputtype = "both", graphname = "2_data/ModelFitControls.pdf")
  save(output2, file="2_data/Posteriors_controlData.RData")
}