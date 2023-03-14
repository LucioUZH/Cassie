# ==============================================================================
# Cassandra Gunasekaram, AIM UZH, PhD - Chimp genetics project, 21.02.2023
# Analyse model/test outputs
# Last update: 14.03.2023
# ==============================================================================

# ______________________________________________________________________________
# Workspace preparations and data input----

  rm(list=ls())
  lapply(names(sessionInfo()$otherPkgs), function(pkgs) detach(paste0('package:',pkgs),character.only = T,unload = T,force=T))
  setwd("C:/Users/cassi/OneDrive - Universität Zürich UZH/PhD/Genetics/Analysis/Repository")
  library(dplyr)
  library(ggpubr)

## GENERAL FUNCTIONS............................................................
# Not in function
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
# ----
# ______________________________________________________________________________
# Bayesian models-----
  
  library(brms)
  library(rstan)
  library(tidybayes)
  library(posterior)
  library(ggforce)
  
  load("Results/bayesiandefault.RData")
  #load("Results/bayesianwide.RData")

# ANALYSE MODEL OUTPUTS........................................................
  
# View model outputs
  summary(bm.complex.maxibd)

  draws <- as_draws_array(bm.complex.snp)
  summarise_draws(draws, default_summary_measures())
  
  get_variables(bm.complex.snp)
  bm.complex.snp %>%
    spread_draws(b_Intercept, b_z.log.similarity) %>%
    mutate(snp_mean = b_Intercept + b_z.log.similarity) %>%
  ggplot(aes(snp_mean, fill = after_stat(x > 0))) +
    stat_halfeye() +
    #geom_vline(xintercept = c(-.8, .8), linetype = "dashed") +
    scale_fill_manual(values = c("gray80", "skyblue"))
  
  posterior_summary(bm.complex.snp)

  #plot(bm.complex.snp, which = 1:2)
  
  mcmc_plot(bm.complex.snp, type = "acf_bar")

## Figure XY....................................................................
  
# Function to extract CI's
  ci <- function(model, cred){
    output <- bayestestR::ci(
      model,
      ci = cred,
      method = "ETI",
      effects = c("fixed", "random", "all"),
      component = c("conditional", "zi", "zero_inflated", "all"),
      parameters = NULL,
      verbose = TRUE,
      BF = 1)
    return(output)
  }
  
# Function to store model outputs in graphdata
  compilegraph <- function(data, outcome, predictor, model){
    data[nrow(data)+1,] <- c(outcome, predictor, fixef(model)[2,1],
                             paste0(outcome,predictor,"67"),ci(model,0.67)[2, 3],
                             "67")
    data[nrow(data)+1,] <- c(outcome, predictor, fixef(model)[2,1],
                             paste0(outcome,predictor,"67"),ci(model,0.67)[2, 4],
                             "67")
    data[nrow(data)+1,] <- c(outcome, predictor, fixef(model)[2,1],
                             paste0(outcome,predictor,"87"),ci(model,0.87)[2, 3],
                             "87")
    data[nrow(data)+1,] <- c(outcome, predictor, fixef(model)[2,1],
                             paste0(outcome,predictor,"87"),ci(model,0.87)[2, 4],
                             "87")
    data[nrow(data)+1,] <- c(outcome, predictor, fixef(model)[2,1],
                             paste0(outcome,predictor,"97"),ci(model,0.97)[2, 3],
                             "97")
    data[nrow(data)+1,] <- c(outcome, predictor, fixef(model)[2,1],
                             paste0(outcome,predictor,"97"),ci(model,0.97)[2, 4],
                             "97")
    return(data)
  }
  
# Create data frame to store
  graphdata <- data.frame(model = character(), predictor = character(), 
                          estimate = numeric(), cilevel = character(),
                          ci = numeric(), nrci = numeric())
  
# Add modes tp graphdata (from lowest to highest in the graph)
  graphdata <- compilegraph(graphdata, "complex", "IBD max", bm.complex.ibdmax)
  graphdata <- compilegraph(graphdata, "complex", "IBD mean", bm.complex.ibdmean)
  graphdata <- compilegraph(graphdata, "complex", "IBD yn", bm.complex.ibdyn)
  graphdata <- compilegraph(graphdata, "complex", "rare alleles", bm.complex.snp)
  
  graphdata <- compilegraph(graphdata, "simple", "IBD max", bm.simple.ibdmax)
  graphdata <- compilegraph(graphdata, "simple", "IBD mean", bm.simple.ibdmean)
  graphdata <- compilegraph(graphdata, "simple", "IBD yn", bm.simple.ibdyn)
  graphdata <- compilegraph(graphdata, "simple", "rare alleles", bm.simple.snp)
  
  graphdata <- compilegraph(graphdata, "tool", "IBD mean", bm.tool.ibdmean)
  graphdata <- compilegraph(graphdata, "tool", "IBD yn", bm.tool.ibdyn)
  graphdata <- compilegraph(graphdata, "tool", "rare alleles", bm.tool.snp)
  
  graphdata <- compilegraph(graphdata, "nontool", "IBD mean", bm.ntool.ibdmean)
  graphdata <- compilegraph(graphdata, "nontool", "IBD yn", bm.ntool.ibdyn)
  graphdata <- compilegraph(graphdata, "nontool", "rare alleles", bm.ntool.snp)
  
  # Add position on y axis and prepare parameters
  ypos <- c(1,2,3,4,7,8,9,10,13,14,15,18,19,20)
  ybreaks <- c(3,9,14,19)
  ypos <- ypos[rep(1:length(ypos),each=6)]
  graphdata$ypos <- ypos
  graphdata$model <- as.factor(graphdata$model)
  graphdata$ci <- as.numeric(graphdata$ci)
  graphdata$estimate <- as.numeric(graphdata$estimate)
  graphdata$predictor <- factor(graphdata$predictor,
                                levels = c("rare alleles", "IBD yn", "IBD mean", "IBD max"))
  levels(graphdata$predictor) <- c("Rare alleles", "IBD yn", "Mean IBD length", "Max IBD length")
  
  
  # Graph
  ggplot(data = graphdata)  +
    geom_link2(aes(x = ci, y = ypos, col = predictor, group = cilevel,
                   size = nrci, alpha = nrci), lineend = "round")+
    geom_point(aes(x = estimate, y = ypos, col = predictor), size = 2)+ # dot for legend
    geom_point(aes(x = estimate, y= ypos, fill = predictor), size = 4,
               shape = 21, col = "black", show.legend = F)+
    labs(size = "CI", col = "Genetic predictor")+
    scale_fill_manual(values = c("orchid1", "deepskyblue", "aquamarine2", "goldenrod"))+
    scale_size_manual(values = c(2.5,1.5,0.9))+
    scale_alpha_manual(values = c(1,0.7,0.5), guide = "none") +
    scale_colour_manual(values = c("orchid1", "deepskyblue", "aquamarine2", "goldenrod"))+
    
    geom_vline(xintercept =  0, lty = 2, size = 0.5, col = "black") +
    geom_hline(yintercept =  9.5, lty = 2, size = 0.3, col = "grey") +
    
    xlab("Estimate ± CI [97%, 87%, 67%]")+
    xlim(-2.2,2.2)+
    scale_y_continuous(breaks = ybreaks, label = c(
      "Composite", "Simple", 
      "Tool use","Nontool"))+
    labs(size = "CI", col = "Genetic predictor")+
    
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white", 
                                          colour = "black", linetype = "solid"),
          panel.grid.major.x = element_line(size = 0.2, color = "grey",
                                            linetype = 2),
          panel.grid.minor.x = element_line(size = 0.2, color = "grey",
                                            linetype = 2),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.key = element_blank())
  
  
## CLEAR WORKSPACE..............................................................s
  rm(list= ls()[!(ls() == "%!in%")])  
  detach("package:brms", unload=TRUE)
  detach("package:rstan", unload=TRUE)
  detach("package:tidybayes", unload=TRUE)
  detach("package:posterior", unload=TRUE)
  detach("package:ggforce", unload=TRUE)
  
# ----
# ______________________________________________________________________________
# Mantel correlation tests & correlogram----
  
  library(corrplot)
  
  load("Results/manteltests.RData")
  
## VISUALISE SIMPLE MANTEL TEST.................................................
  corrplot(mantel.all.r, p.mat = mantel.all.p, method = 'pie', 
           type = 'lower', insig='p-value',sig.level = -1, 
           diag=FALSE, tl.col = 1) 
  
## VISUALISE PARTIAL MANTEL TESTS...............................................
  
  outer = FALSE
  line = -2
  cex = 1.2
  adj  = 0.025
  
  par(mfrow=c(1,2), font.main=2)
  corrplot(mantel.snp.r, p.mat = mantel.snp.p, method = 'pie', type = 'upper', 
           insig='p-value', sig.level = -1, diag=FALSE, tl.col = 1)
  title(outer=outer,adj=adj,main="A) Controlling for Rare alleles",cex.main=cex,col="black",font=2,line=line)
  corrplot(mantel.ibd.r, p.mat = mantel.ibd.p, method = 'pie', type = 'upper', 
           insig='p-value', sig.level = -1, diag=FALSE, tl.col = 1)
  title(outer=outer,adj=adj,main="B) Controlling for IBD",cex.main=cex,col="black",font=2,line=line)
  
## VISUALISE SIMPLE MANTEL TESTS OF COMPLEX TOOLS...............................
  par(mfrow=c(2,2), font.main=2)
  corrplot(mantel.accstone.r, p.mat = mantel.accstone.p, method = 'pie', type = 'upper', 
           insig='p-value', sig.level = -1, diag=FALSE, tl.col = 1)
  title(outer=outer,adj=adj,main="A) Stone Throwing",cex.main=cex,col="black",font=2,line=line)
  corrplot(mantel.honey.r, p.mat = mantel.honey.p, method = 'pie', type = 'upper', 
           insig='p-value', sig.level = -1, diag=FALSE, tl.col = 1)
  title(outer=outer,adj=adj,main="B) Honey",cex.main=cex,col="black",font=2,line=line)
  corrplot(mantel.nut.r, p.mat = mantel.nut.p, method = 'pie', type = 'upper', 
           insig='p-value', sig.level = -1, diag=FALSE, tl.col = 1)
  title(outer=outer,adj=adj,main="C) Nut cracking",cex.main=cex,col="black",font=2,line=line)
  corrplot(mantel.underground.r, p.mat = mantel.underground.p, method = 'pie', type = 'upper', 
           insig='p-value', sig.level = -1, diag=FALSE, tl.col = 1)
  title(outer=outer,adj=adj,main="D) Underground",cex.main=cex,col="black",font=2,line=line)
  

# Simple Mantel simple tools
  par(mfrow=c(2,2), font.main=2)
  corrplot(mantel.comm.r, p.mat = mantel.comm.p, method = 'pie', type = 'upper', 
           insig='p-value', sig.level = -1, diag=FALSE, tl.col = 1)
  title(outer=outer,adj=adj,main="A) Communication",cex.main=cex,col="black",font=2,line=line)
  corrplot(mantel.insect.r, p.mat = mantel.insect.p, method = 'pie', type = 'upper', 
           insig='p-value', sig.level = -1, diag=FALSE, tl.col = 1)
  title(outer=outer,adj=adj,main="B) Insect foraging",cex.main=cex,col="black",font=2,line=line)
  corrplot(mantel.foraging.r, p.mat = mantel.foraging.p, method = 'pie', type = 'upper', 
           insig='p-value', sig.level = -1, diag=FALSE, tl.col = 1)
  title(outer=outer,adj=adj,main="C) Foraging",cex.main=cex,col="black",font=2,line=line)
  
## CORRELOGRAMS FOR SPATIAL AUTOCORRELATION.....................................
  
  library(ecodist)
  
  load("Results/similaritymatrices.RData")
  
# Complexity groups
  par(mfrow=c(2,2), font.main=4, mar=c(5,4,2,2))
  plot(mgram(cultdist0, geodist), 
       main='(a) Nontool', xlab='Distance / km')
  plot(mgram(cultdist1, geodist), 
       main='(b) Simple', xlab='Distance / km')
  plot(mgram(cultdist2, geodist), 
       main='(c) Complex', xlab='Distance / km')
  
# Complex tool groups
  par(mfrow=c(2,2), font.main=4, mar=c(5,4,2,2))
  plot(mgram(cultdistc, geodist), 
       main='(A) Acc. stone throwing', xlab='Distance / km')
  plot(mgram(cultdisth, geodist), 
       main='(b) Honey toolset', xlab='Distance / km')
  plot(mgram(cultdistn, geodist), 
       main='(c) Nut cracking', xlab='Distance / km')
  plot(mgram(cultdistu, geodist), 
       main='(D) Underground', xlab='Distance / km')
  
## CLEAR WORKSPACE..............................................................s
  rm(list= ls()[!(ls() %in% c('ddyads.raw','ddyads', "%!in%"))])

# ---- 

