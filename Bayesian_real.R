################################################################################
# CASSANDRA GUNASEKARAM, BAYESIAN REGRESSIONS MAIN #############################
# 15.07.2024 ###################################################################

library(dplyr)
library(cgwtools)

rm(list=ls())

setwd("C:/Users/cassi/OneDrive - Universität Zürich UZH/PhD/Project 1 Genes and Culture in Chimpanzees/Sorted")

##################################################################     LOAD DATA -----

  # IMPORT DATA
    d.dyads <- read.csv("Data/Processed/d.dyads.csv", sep = ";")
  
  # PREPARE VARIABLES
    str(d.dyads)
    d.dyads[d.dyads == "na"] <- NA
    d.dyads <- filter(d.dyads, !is.na(shared))
    
    d.dyads$z.log.nepra <- scale(log(d.dyads$neprasimilarity))
    d.dyads$z.obstime <- scale(d.dyads$min_observationtime_months)
    d.dyads$behaviour <- as.factor(d.dyads$behaviour)
    
################################################################################
## MODEL FITTING
################################################################################
########################################     FUNCTIONS AND PRESETS TO RUN MODELS -----

  # LIBRARIES
    library(rstan)
    library(brms)
    
  # SHORTCUT TO CHANGE ALL AT THE SAME TIME
    # ctrl+Alt+shift+M
  
  # !IN 
    '%!in%' <- function(x,y)!('%in%'(x,y))
  
  # PRIORS (WEAK)
    priors <- c(prior(normal(0,1), class=Intercept),
                prior(normal(0,1), class=b),
                prior(exponential(1), class=sd))
  
  # FUNCTION TO RUN BAYESIAN MODEL WITH 2,000 ITERATIONS AND 1,000 WARMUP
    get_bayesianmodel <- function(form, fam, dataset, priors){
      bm <- brm(formula = form,family = fam,data = dataset,
                prior= priors, 
                cores = 16, seed = 100, 
                warmup=1000, iter=2000, chains=4,
                control=list(adapt_delta =0.99))}
    
    
################################################################### NEPRA MODELS -----
    
  # SUBSET DATA TO NON-TOOL, SIMPLE AND COMPLEX
    d.nontool <- droplevels(filter(d.dyads, complexity == 0))
    d.simple <- droplevels(filter(d.dyads, complexity == 1))
    d.complex <- droplevels(filter(d.dyads, complexity == 2))
    
  # FORMULAS WITH HABITAT
    fneprab <- bf(shared ~ neprabinary + samehabitat + z.obstime + (1|dyadID) + (1|behaviour))
    fneprasim <- bf(shared ~ z.log.nepra + samehabitat + z.obstime + (1|dyadID) + (1|behaviour))
    
  # NEPRA BINARY MODELS
    m.c.neprab.hab <- get_bayesianmodel(fneprab, "bernoulli", d.complex, priors)
    m.s.neprab.hab <- get_bayesianmodel(fneprab, "bernoulli", d.simple, priors)
    m.n.neprab.hab <- get_bayesianmodel(fneprab, "bernoulli", d.nontool, priors)
    
  # NEPRA SIMILARITY MODELS
    m.c.neprasim.hab <- get_bayesianmodel(fneprasim, "bernoulli", d.complex, priors)
    m.s.neprasim.hab <- get_bayesianmodel(fneprasim, "bernoulli", d.simple, priors)
    m.n.neprasim.hab <- get_bayesianmodel(fneprasim, "bernoulli", d.nontool, priors)
    
  # FORMULAS WITH DISTANCE TO REFUGIUM
    fneprab <- bf(shared ~ neprabinary + z.sqrt.ref + z.obstime + (1|dyadID) + (1|behaviour))
    fneprasim <- bf(shared ~ z.log.nepra + z.sqrt.ref + z.obstime + (1|dyadID) + (1|behaviour))
    
  # NEPRA BINARY MODELS
    m.c.neprab.ref <- get_bayesianmodel(fneprab, "bernoulli", d.complex, priors)
    m.s.neprab.ref <- get_bayesianmodel(fneprab, "bernoulli", d.simple, priors)
    m.n.neprab.ref <- get_bayesianmodel(fneprab, "bernoulli", d.nontool, priors)
    
  # NEPRA SIMILARITY MODELS
    m.c.neprasim.ref <- get_bayesianmodel(fneprasim, "bernoulli", d.complex, priors)
    m.s.neprasim.ref <- get_bayesianmodel(fneprasim, "bernoulli", d.simple, priors)
    m.n.neprasim.ref <- get_bayesianmodel(fneprasim, "bernoulli", d.nontool, priors)
    
  
  # FORMULAS WITH ANNUAL PRECIPITATION
    fneprab <- bf(shared ~ neprabinary + z.log.precip + z.obstime + (1|dyadID) + (1|behaviour))
    fneprasim <- bf(shared ~ z.log.nepra + z.log.precip + z.obstime + (1|dyadID) + (1|behaviour))
    
  # NEPRA BINARY MODELS
    m.c.neprab.prec <- get_bayesianmodel(fneprab, "bernoulli", d.complex, priors)
    m.s.neprab.prec <- get_bayesianmodel(fneprab, "bernoulli", d.simple, priors)
    m.n.neprab.prec <- get_bayesianmodel(fneprab, "bernoulli", d.nontool, priors)
    
  # NEPRA SIMILARITY MODELS
    m.c.neprasim.prec <- get_bayesianmodel(fneprasim, "bernoulli", d.complex, priors)
    m.s.neprasim.prec <- get_bayesianmodel(fneprasim, "bernoulli", d.simple, priors)
    m.n.neprasim.prec <- get_bayesianmodel(fneprasim, "bernoulli", d.nontool, priors)
    
    
##################################################################### IBD MODELS -----
    
  # SUBSET DATA TO NON-TOOL, SIMPLE AND COMPLEX EXCLUDING OUTAMBA KILIMI (561 DYADS)
    d.nontool <- droplevels(filter(d.dyads, site1 != "Outamba-Kilimi" & site2!= "Outamba-Kilimi" & complexity == 0))
    d.simple <- droplevels(filter(d.dyads, site1 != "Outamba-Kilimi" & site2!= "Outamba-Kilimi" & complexity == 1))
    d.complex <- droplevels(filter(d.dyads, site1 != "Outamba-Kilimi" & site2!= "Outamba-Kilimi" & complexity == 2))
    
  # FORMULA FOR HABITAT
    fibd <- bf(shared ~ ibdbinary + samehabitat + z.obstime + (1|dyadID) + (1|behaviour))
    
  # IBD MODELS
    m.c.ibd.hab <- get_bayesianmodel(fibd, "bernoulli", d.complex, priors)
    m.s.ibd.hab <- get_bayesianmodel(fibd, "bernoulli", d.simple, priors)
    m.n.ibd.hab <- get_bayesianmodel(fibd, "bernoulli", d.nontool, priors)

  # FORMULA FOR DISTANCE TO REFUGIUM
    fibd <- bf(shared ~ ibdbinary + z.sqrt.ref + z.obstime + (1|dyadID) + (1|behaviour))
    
  # IBD MODELS
    m.c.ibd.ref <- get_bayesianmodel(fibd, "bernoulli", d.complex, priors)
    m.s.ibd.ref <- get_bayesianmodel(fibd, "bernoulli", d.simple, priors)
    m.n.ibd <- get_bayesianmodel(fibd, "bernoulli", d.nontool, priors)
    
  # FORMULA FOR ANNUAL PRECIPITATION
    fibd <- bf(shared ~ ibdbinary + z.log.precip + z.obstime + (1|dyadID) + (1|behaviour))
    
  # IBD MODELS
    m.c.ibd.prec <- get_bayesianmodel(fibd, "bernoulli", d.complex, priors)
    m.s.ibd.prec <- get_bayesianmodel(fibd, "bernoulli", d.simple, priors)
    m.n.ibd.prec <- get_bayesianmodel(fibd, "bernoulli", d.nontool, priors)
    
#################################################################### SAVE MODELS -----
    
    save(m.c.neprab.hab, m.c.neprasim.hab, m.c.ibd.hab,
         m.c.neprab.hab, m.c.neprasim.hab, m.c.ibd.hab,
         m.c.neprab.hab, m.c.neprasim.hab, m.c.ibd.hab,
         file = "Results/Bayesian/models_habitat.RData")
    save(m.c.neprab.ref, m.c.neprasim.ref, m.c.ibd.ref,
         m.c.neprab.ref, m.c.neprasim.ref, m.c.ibd.ref,
         m.c.neprab.ref, m.c.neprasim.ref, m.c.ibd.ref,
         file = "Results/Bayesian/models_refugium.RData")
    save(m.c.neprab.prec, m.c.neprasim.prec, m.c.ibd.prec,
         m.c.neprab.prec, m.c.neprasim.prec, m.c.ibd.prec,
         m.c.neprab.precb, m.c.neprasim.prec, m.c.ibd.prec,
         file = "Results/Bayesian/models_precipitation.RData")
    
################################################################################
## MODEL DIAGNOSTICS
################################################################################
################################################################ MODEL ESTIMATES ------
  
  # LIBRARIES
    library(grid)
    library(gridExtra)
    
  # FUNCTION TO GET MODEL SUMMARY
    get_summary <- function(model_name, model){
      modelfit <- as.data.frame(summary(model$fit))
      posterior_draws <- as.data.frame(model$fit)
      df <- data.frame(modelname = model_name,
                       variables = rownames(modelfit)[1:6],
                       estimate.mean = modelfit[1:6,1], 
                       estimate.sd = modelfit[1:6,3], 
                       estimate.ci2.5 = modelfit[1:6,4],
                       estimate.ci97.5 = modelfit[1:6,8], 
                       Post0 = NA, 
                       OR = exp(modelfit[1:6,1]), 
                       OR2.5 = exp(modelfit[1:6,4]), 
                       OR97.5 = exp(modelfit[1:6,8]),
                       ESS = modelfit[1:6,9],
                       rhat = modelfit[1:6,10])
      for ( i in 1:6){
        post_draws <- posterior_draws[, i]
        df$Post0[i] <- mean(post_draws > 0)
        
      }
      return(df)
    }
    
  # GET SUMMARY FOR ALL MODELS HABITAT AND PRINT TO PDF
    pdf("Results/Bayesian/md_summaries_habitat.pdf", width = 17, height = 7.6)
    
    
    grid.table(rbind(get_summary("Complex_NePRA_binary", m.c.neprab.hab),
                     get_summary("Simple_NePRA_binary", m.s.neprab.hab),
                     get_summary("Nontool_NePRA_binary", m.n.neprab.hab)))
    
    grid.newpage()
    grid.table(rbind(get_summary("Complex_NePRA_similarity", m.c.neprasim.hab),
                     get_summary("Simple_NePRA_similarity", m.s.neprasim.hab),
                     get_summary("Nontool_NePRA_similarity", m.n.neprasim.hab)))
    
    grid.newpage()
    grid.table(rbind(get_summary("Complex_IBD", m.c.ibd.hab),
                     get_summary("Simple_IBD", m.s.ibd.hab),
                     get_summary("Nontool_IBD", m.n.ibd.hab)))
    
    dev.off()
    
##################################################### POSTERIOR PREDICTIVE CHECK -----
  
  # PREPARATION OF THEME OF PLOTS
    themeset <- theme(
      text = element_text(family = "sans"),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = c(0.8,0.8),
      axis.ticks.y = element_blank(),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.2, 'cm')
    )
  
  # FUNCTION TO GET PLOT
  
    get_ppcheck <- function(title, model){
      pp_check(model, type = "bars") + themeset + 
        scale_x_continuous(breaks = c(0,1), labels = c("0", "1")) + 
        ggtitle(title)
    }
    
    pp.check.hab <- 
      ggarrange(get_ppcheck("Complex ~ NePRA binary", m.c.neprab.hab),
                get_ppcheck("Simple ~ NePRA binary", m.c.neprab.hab),
                get_ppcheck("Non-tool ~ NePRA binary", m.c.neprab.hab),
                get_ppcheck("Complex ~ NePRA similarity", m.c.neprasim.hab),
                get_ppcheck("Simple ~ NePRA similarity", m.c.neprasim.hab),
                get_ppcheck("Non-tool ~ NePRA similarity", m.c.neprasim.hab),
                get_ppcheck("Complex ~ IBD binary", m.c.ibd.hab),
                get_ppcheck("Simple ~ IBD binary", m.c.ibd.hab),
                get_ppcheck("Non-tool ~ IBD binary", m.c.ibd.hab),
                nrow = 3, ncol = 3)

    ggsave(plot = pp.check.hab, file = "Results/Bayesian/md_ppcheck_habitat.jpeg", 
           dpi = 300, width = 15, height = 15, unit = "cm")   

  
########################################################## DENSITY & TRACE PLOTS  -----
    
  # LIBRARIES
    library(bayesplot)
    library(patchwork)
    
  # FUNCTION TO GET PLOTS
    get_plots <- function(title, model, trace_var){
      p.trace <- mcmc_trace(model, pars = trace_var, 
                            facet_args = list(ncol = 1, strip.position = "left")) + 
      labs(
        title = " ",
        subtitle = "Model convergence"
      )
      
      # Create marginal posterior distributions plot for fixed and randomeffects
      p.posterior <- mcmc_areas(model, pars = trace_var) +
        labs(
          title = title,
          subtitle = "Posterior distribution with medians and 95% intervals"
        )
    
      # Combine all plots using patchwork
      combined_plot <- p.posterior |p.trace
      return(combined_plot) 
    }
    
    
  # GET PLOTS AND SAVE AS PDF
    pdf("Results/Bayesian/md_tracedensity_habitat.pdf", width = 15, height = 10)
    

    # Arrange pages
    get_plots("Complex ~ NePRA binary", m.c.neprab.hab,
              c("b_Intercept", "b_neprabinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Simple ~ NePRA binary", m.s.neprab.hab,
              c("b_Intercept", "b_neprabinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Non-tool ~ NePRA binary", m.n.neprab.hab,
              c("b_Intercept", "b_neprabinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Complex ~ NePRA similarity", m.c.neprasim.hab,
              c("b_Intercept", "b_z.log.nepra", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Simple ~ NePRA similarity", m.s.neprasim.hab,
              c("b_Intercept", "b_z.log.nepra", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept")) 
    get_plots("Non-tool ~ NePRA similarity", m.n.neprasim.hab,
              c("b_Intercept", "b_z.log.nepra", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Complex ~ IBD binary", m.c.ibd.hab,
              c("b_Intercept", "b_ibdbinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Simple ~ IBD binary", m.s.ibd.hab,
              c("b_Intercept", "b_ibdbinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept")) 
    get_plots("Non-tool ~ IBD binary", m.n.ibd.hab,
              c("b_Intercept", "b_ibdbinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept")) 
    dev.off()
    

######################################################### TEST OUTLIERS WITH LOO -----
    
    
  # FUNCTION TO GET LOO PLOT
    get_pareto <- function(title, model){
      loo1 <- loo(model)
      pareto1 <- loo1$diagnostics$pareto_k
      plot(pareto1, type = "p", ylim = c(0, 0.8),
           xlab = "Data Point Index", ylab = "Pareto k Value",
           main = title)
      abline(h = 0.7, lty = 2, col = "red")
    }

  # GET PLOTS
    
    pdf("Results/Bayesian/md_loo_habitat.pdf", width = 8, height = 7)
    
    par(mfrow=c(3,3))
    
    get_pareto("Complex ~ NePRA binary", m.c.neprab.hab)
    get_pareto("Complex ~ NePRA similarity", m.c.neprasim.hab)
    get_pareto("Complex ~ IBD binary", m.c.ibd.hab)
    get_pareto("Simple ~ NePRA binary", m.s.neprab.hab)
    get_pareto("Simple ~ NePRA similarity", m.s.neprasim.hab)
    get_pareto("Simple ~ IBD binary", m.s.ibd.hab)
    get_pareto("Non-tool ~ NePRA binary", m.n.neprab.hab)
    get_pareto("Non-tool ~ NePRA similarity", m.n.neprasim.hab)
    get_pareto("Non-tool ~ IBD binary", m.n.ibd.hab)
    
    dev.off()
      
  # CHECK MODELS THAT HAVE PARETO > 0.7
    
  # 1. Compex ~ IBD
    d.complex <- droplevels(filter(d.dyads, site1 != "Outamba-Kilimi" & site2!= "Outamba-Kilimi" & complexity == 2))
    d.complex2 <- d.complex[-which(loo(m.c.ibd.hab)$diagnostics$pareto_k>0.7),]

    # Run model without these data points:
    m.c.ibd.LOO <- get_bayesianmodel(fibd, "bernoulli", d.complex2, priors)
    summary(m.c.ibd.hab)
    summary(m.c.ibd.LOO)
    
  # 2. Nontool ~ NePRAb
    d.nontool <- droplevels(filter(d.dyads, complexity == 0))
    d.nontool2 <- d.nontool[-which(loo(m.n.neprab.hab)$diagnostics$pareto_k>0.7),]
    
    # Run model without these data points:
    m.n.neprab.LOO <- get_bayesianmodel(fneprab, "bernoulli", d.nontool2, priors)
    summary(m.n.neprab.hab)
    summary(m.n.neprab.LOO)
    
  # 3. Nontool ~ NePRAsim
    d.nontool <- droplevels(filter(d.dyads, complexity == 0))
    d.nontool2 <- d.nontool[-which(loo(m.n.neprasim.hab)$diagnostics$pareto_k>0.7),]
    
    # Run model without these data points:
    m.n.neprasim.LOO <- get_bayesianmodel(fneprasim, "bernoulli", d.nontool2, priors)
    summary(m.n.neprasim.hab)
    summary(m.n.neprasim.LOO)

  # 4. Non-tool ~ IBD
    d.nontool <- droplevels(filter(d.dyads, site1 != "Outamba-Kilimi" & site2!= "Outamba-Kilimi" & complexity == 0))
    d.nontool2 <- d.nontool[-which(loo(m.n.ibd.hab)$diagnostics$pareto_k>0.7),]
    
    # Run model without these data points:
    m.n.ibd.LOO <- get_bayesianmodel(fibd, "bernoulli", d.nontool2, priors)
    summary(m.n.ibd.hab)
    summary(m.n.ibd.LOO)
