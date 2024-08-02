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
    get_bayesian_model <- function(formel, df.culture, df.corrmatrix, varlist) {
      
      # Bayesian model
      bm <- brm(formula = formel,
               family = "bernoulli", data = df.culture, data2 = list(A = df.corrmatrix), prior = priors,
               cores = 8, seed = 100, warmup=1000, iter=2000, chains=4, control=list(adapt_delta =0.99))
    }
  

#################### MATRIX WITH CORRELATION COEFFICIENT OF 0.5 FOR NEPRA MODELS -----

  # FUNCTION TO GET POPULATIONS OF EACH DYAD
    get_pops <- function(dyad) {
      strsplit(dyad, "_")[[1]]
    }
  
  # GET ALL DYADS
    l.dyads <- unique(sort(d.dyads$dyadID))
    n <- length(l.dyads)
  
  # CREATE MATRIX
    d.covmatrix.empty <- matrix(0, nrow = n, ncol = n)
    
    rownames(d.covmatrix.empty) <- l.dyads
    colnames(d.covmatrix.empty) <- l.dyads
  
  # CORRELATION COEFFICIENT
    summ <- 0.5
  
  # FILL MATRIX
    d.covmatrix <- d.covmatrix.empty
    for (i in 1:n) {
      for (j in 1:n) {
        # Get the names from the row and column
        popij <- get_pops(rownames(d.covmatrix)[i])
        popkl <- get_pops(colnames(d.covmatrix)[j])
        
        # Check if both first and second names match
        if (popij[1] == popkl[1] && popij[2] == popkl[2]) {
          d.covmatrix[i, j] <- d.covmatrix[i, j] + 1.000000002
        } else if (popij[1] == popkl[1] || popij[2] == popkl[2]) {
          d.covmatrix[i, j] <- d.covmatrix[i, j] + summ
        } else if (popij[1] == popkl[2] || popij[2] == popkl[1]) {
          d.covmatrix[i, j] <- d.covmatrix[i, j] - summ
        }
      }
    }
  
  # REPLACE ALL NEGATIVE VALUES BY THE POSITIVE
    d.covmatrix[d.covmatrix == -summ] <- summ
  
  
  # VIEW MATRIX IN HEATMAP
    # Convert to long format
    # melted_data <- reshape2::melt(d.covmatrix)
    # table(melted_data$value)
    # my_colors<- colorRampPalette(c("cyan", "white", "deeppink3"))
    # heatmap(d.covmatrix, Rowv = NA, Colv = NA, col = my_colors(4))
    # legend(x="right", legend=c( "", "-0.5", "0", "0.5", "1"),fill=my_colors(4))

################################################################### NEPRA MODELS -----

  # SUBSET DATA TO NON-TOOL, SIMPLE AND COMPLEX
    d.nontool <- droplevels(filter(d.dyads, complexity == 0))
    d.simple <- droplevels(filter(d.dyads, complexity == 1))
    d.complex <- droplevels(filter(d.dyads, complexity == 2))
  
  # FORMULAS WITH HABITAT
    fneprab <- bf(shared ~ neprabinary + samehabitat + z.obstime + (1 | gr(dyadID, cov = A)) + (1 | behaviour))
    fneprasim <- bf(shared ~ z.log.nepra + samehabitat + z.obstime + (1 | gr(dyadID, cov = A)) + (1 | behaviour))
  
  # NEPRA BINARY MODELS
    m.c.neprab.hab.cov0.5 <- get_bayesianmodel(fneprab, "bernoulli", d.complex, priors)
    m.s.neprab.hab.cov0.5 <- get_bayesianmodel(fneprab, "bernoulli", d.simple, priors)
    m.n.neprab.hab.cov0.5 <- get_bayesianmodel(fneprab, "bernoulli", d.nontool, priors)

  # NEPRA SIMILARITY MODELS
    m.c.neprasim.hab.cov0.5 <- get_bayesianmodel(fneprasim, "bernoulli", d.complex, priors)
    m.s.neprasim.hab.cov0.5 <- get_bayesianmodel(fneprasim, "bernoulli", d.simple, priors)
    m.n.neprasim.hab.cov0.5 <- get_bayesianmodel(fneprasim, "bernoulli", d.nontool, priors)

###################### MATRIX WITH CORRELATION COEFFICIENT OF 0.5 FOR IBD MODELS -----
    
  # GET ALL DYADS WITHOUT OUTAMBA-KILIMI
    l.dyads <- unique(sort(d.dyads[d.dyads$site1 != "Outamba-Kilimi" & d.dyads$site2 != "Outamba-Kilimi",]$dyadID))
    n <- length(l.dyads)
    
    # CREATE MATRIX
    d.covmatrix.empty.ibd <- matrix(0, nrow = n, ncol = n)
    
    rownames(d.covmatrix.empty.ibd) <- l.dyads
    colnames(d.covmatrix.empty.ibd) <- l.dyads
    
    # CORRELATION COEFFICIENT
    summ <- 0.5
    
    # FILL MATRIX
    d.covmatrix <- d.covmatrix.empty.ibd
    for (i in 1:n) {
      for (j in 1:n) {
        # Get the names from the row and column
        popij <- get_pops(rownames(d.covmatrix)[i])
        popkl <- get_pops(colnames(d.covmatrix)[j])
        
        # Check if both first and second names match
        if (popij[1] == popkl[1] && popij[2] == popkl[2]) {
          d.covmatrix[i, j] <- d.covmatrix[i, j] + 1.000000002
        } else if (popij[1] == popkl[1] || popij[2] == popkl[2]) {
          d.covmatrix[i, j] <- d.covmatrix[i, j] + summ
        } else if (popij[1] == popkl[2] || popij[2] == popkl[1]) {
          d.covmatrix[i, j] <- d.covmatrix[i, j] - summ
        }
      }
    }
    
    # REPLACE ALL NEGATIVE VALUES BY THE POSITIVE
    d.covmatrix[d.covmatrix == -summ] <- summ
    
##################################################################### IBD MODELS -----
    
  # SUBSET DATA TO NON-TOOL, SIMPLE AND COMPLEX
    d.nontool <- droplevels(filter(d.dyads, site1 != "Outamba-Kilimi" & site2!= "Outamba-Kilimi" & complexity == 0))
    d.simple <- droplevels(filter(d.dyads, site1 != "Outamba-Kilimi" & site2!= "Outamba-Kilimi" & complexity == 1))
    d.complex <- droplevels(filter(d.dyads, site1 != "Outamba-Kilimi" & site2!= "Outamba-Kilimi" & complexity == 2))
    
    
    # FORMULA WITH HABITAT
    fibd <- bf(shared ~ ibdbinary + samehabitat + z.obstime + (1 | gr(dyadID, cov = A)) + (1 | behaviour))

    # NEPRA BINARY MODELS
    m.c.ibd.hab.cov0.5 <- get_bayesianmodel(fibd, "bernoulli", d.complex, priors)
    m.s.ibd.hab.cov0.5 <- get_bayesianmodel(fibd, "bernoulli", d.simple, priors)
    m.n.ibd.hab.cov0.5 <- get_bayesianmodel(fibd, "bernoulli", d.nontool, priors)
    

#################### MATRIX WITH CORRELATION COEFFICIENT OF 0.1 FOR NEPRA MODELS -----
    
  # CORRELATION COEFFICIENT OF 0.1
    summ <- 0.1
    
  # FILL MATRIX
    d.covmatrix <- d.covmatrix.empty
    for (i in 1:n) {
      for (j in 1:n) {
        # Get the names from the row and column
        popij <- get_pops(rownames(d.covmatrix)[i])
        popkl <- get_pops(colnames(d.covmatrix)[j])
        
        # Check if both first and second names match
        if (popij[1] == popkl[1] && popij[2] == popkl[2]) {
          d.covmatrix[i, j] <- d.covmatrix[i, j] + 1.000000002
        } else if (popij[1] == popkl[1] || popij[2] == popkl[2]) {
          d.covmatrix[i, j] <- d.covmatrix[i, j] + summ
        } else if (popij[1] == popkl[2] || popij[2] == popkl[1]) {
          d.covmatrix[i, j] <- d.covmatrix[i, j] - summ
        }
      }
    }
    
  # REPLACE ALL NEGATIVE VALUES BY THE POSITIVE
    d.covmatrix[d.covmatrix == -summ] <- summ
    
################################################################### NEPRA MODELS -----
    
  # SUBSET DATA TO NON-TOOL, SIMPLE AND COMPLEX
    d.nontool <- droplevels(filter(d.dyads, complexity == 0))
    d.simple <- droplevels(filter(d.dyads, complexity == 1))
    d.complex <- droplevels(filter(d.dyads, complexity == 2))
    
  # FORMULAS WITH HABITAT
    fneprab <- bf(shared ~ neprabinary + samehabitat + z.obstime + (1 | gr(dyadID, cov = A)) + (1 | behaviour))
    fneprasim <- bf(shared ~ z.log.nepra + samehabitat + z.obstime + (1 | gr(dyadID, cov = A)) + (1 | behaviour))
    
  # NEPRA BINARY MODELS
    m.c.neprab.hab.cov0.1 <- get_bayesianmodel(fneprab, "bernoulli", d.complex, priors)
    m.s.neprab.hab.cov0.1 <- get_bayesianmodel(fneprab, "bernoulli", d.simple, priors)
    m.n.neprab.hab.cov0.1 <- get_bayesianmodel(fneprab, "bernoulli", d.nontool, priors)
    
  # NEPRA SIMILARITY MODELS
    m.c.neprasim.hab.cov0.1 <- get_bayesianmodel(fneprasim, "bernoulli", d.complex, priors)
    m.s.neprasim.hab.cov0.1 <- get_bayesianmodel(fneprasim, "bernoulli", d.simple, priors)
    m.n.neprasim.hab.cov0.1 <- get_bayesianmodel(fneprasim, "bernoulli", d.nontool, priors)
    
###################### MATRIX WITH CORRELATION COEFFICIENT OF 0.5 FOR IBD MODELS -----

  # CORRELATION COEFFICIENT
    summ <- 0.1
    
  # FILL MATRIX
    d.covmatrix <- d.covmatrix.empty.ibd
    for (i in 1:n) {
      for (j in 1:n) {
        # Get the names from the row and column
        popij <- get_pops(rownames(d.covmatrix)[i])
        popkl <- get_pops(colnames(d.covmatrix)[j])
        
        # Check if both first and second names match
        if (popij[1] == popkl[1] && popij[2] == popkl[2]) {
          d.covmatrix[i, j] <- d.covmatrix[i, j] + 1.000000002
        } else if (popij[1] == popkl[1] || popij[2] == popkl[2]) {
          d.covmatrix[i, j] <- d.covmatrix[i, j] + summ
        } else if (popij[1] == popkl[2] || popij[2] == popkl[1]) {
          d.covmatrix[i, j] <- d.covmatrix[i, j] - summ
        }
      }
    }
    
  # REPLACE ALL NEGATIVE VALUES BY THE POSITIVE
    d.covmatrix[d.covmatrix == -summ] <- summ
    
##################################################################### IBD MODELS -----
    
  # SUBSET DATA TO NON-TOOL, SIMPLE AND COMPLEX
    d.nontool <- droplevels(filter(d.dyads, site1 != "Outamba-Kilimi" & site2!= "Outamba-Kilimi" & complexity == 0))
    d.simple <- droplevels(filter(d.dyads, site1 != "Outamba-Kilimi" & site2!= "Outamba-Kilimi" & complexity == 1))
    d.complex <- droplevels(filter(d.dyads, site1 != "Outamba-Kilimi" & site2!= "Outamba-Kilimi" & complexity == 2))
    
    
  # FORMULA WITH HABITAT
    fibd <- bf(shared ~ ibdbinary + samehabitat + z.obstime + (1 | gr(dyadID, cov = A)) + (1 | behaviour))
    
  # NEPRA BINARY MODELS
    m.c.ibd.hab.cov0.1 <- get_bayesianmodel(fibd, "bernoulli", d.complex, priors)
    m.s.ibd.hab.cov0.1 <- get_bayesianmodel(fibd, "bernoulli", d.simple, priors)
    m.n.ibd.hab.cov0.1 <- get_bayesianmodel(fibd, "bernoulli", d.nontool, priors)
    
    

#################################################################### SAVE MODELS -----

save(m.c.neprab.hab.cov0.5, m.c.neprasim.hab.cov0.5, m.c.ibd.hab.cov0.5,
     m.c.neprab.hab.cov0.5, m.c.neprasim.hab.cov0.5, m.c.ibd.hab.cov0.5,
     m.c.neprab.hab.cov0.5, m.c.neprasim.hab.cov0.5, m.c.ibd.hab.cov0.5,
     file = "Results/Bayesian/models_habitat_cov0.5.RData")
save(m.c.neprab.hab.cov0.1, m.c.neprasim.hab.cov0.1, m.c.ibd.hab.cov0.1,
     m.c.neprab.hab.cov0.1, m.c.neprasim.hab.cov0.1, m.c.ibd.hab.cov0.1,
     m.c.neprab.hab.cov0.1, m.c.neprasim.hab.cov0.1, m.c.ibd.hab.cov0.1,
     file = "Results/Bayesian/models_habitat_cov0.1.RData")
    
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
    pdf("Results/Bayesian/md_summaries_habitat.cov0.5.pdf", width = 17, height = 7.6)
    
    
    grid.table(rbind(get_summary("Complex_NePRA_binary_cov0.5", m.c.neprab.hab.cov0.5),
                     get_summary("Simple_NePRA_binary_cov0.5", m.s.neprab.hab.cov0.5),
                     get_summary("Nontool_NePRA_binary_cov0.5", m.n.neprab.hab.cov0.5)))
    
    grid.newpage()
    grid.table(rbind(get_summary("Complex_NePRA_similarity_cov0.5", m.c.neprasim.hab.cov0.5),
                     get_summary("Simple_NePRA_similarity_cov0.5", m.s.neprasim.hab.cov0.5),
                     get_summary("Nontool_NePRA_similarity_cov0.5", m.n.neprasim.hab.cov0.5)))
    
    grid.newpage()
    grid.table(rbind(get_summary("Complex_IBD_cov0.5", m.c.ibd.hab.cov0.5),
                     get_summary("Simple_IBD_cov0.5", m.s.ibd.hab.cov0.5),
                     get_summary("Nontool_IBD_cov0.5", m.n.ibd.hab.cov0.5)))
    
    dev.off()
    
    # GET SUMMARY FOR ALL MODELS HABITAT AND PRINT TO PDF
    pdf("Results/Bayesian/md_summaries_habitat.cov0.1.pdf", width = 17, height = 7.6)
    
    
    grid.table(rbind(get_summary("Complex_NePRA_binary_cov0.1", m.c.neprab.hab.cov0.1),
                     get_summary("Simple_NePRA_binary_cov0.1", m.s.neprab.hab.cov0.1),
                     get_summary("Nontool_NePRA_binary_cov0.1", m.n.neprab.hab.cov0.1)))
    
    grid.newpage()
    grid.table(rbind(get_summary("Complex_NePRA_similarity_cov0.1", m.c.neprasim.hab.cov0.1),
                     get_summary("Simple_NePRA_similarity_cov0.1", m.s.neprasim.hab.cov0.1),
                     get_summary("Nontool_NePRA_similarity_cov0.1", m.n.neprasim.hab.cov0.1)))
    
    grid.newpage()
    grid.table(rbind(get_summary("Complex_IBD_cov0.1", m.c.ibd.hab.cov0.1),
                     get_summary("Simple_IBD_cov0.1", m.s.ibd.hab.cov0.1),
                     get_summary("Nontool_IBD_cov0.1", m.n.ibd.hab.cov0.1)))
    
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
  
  # MODELS WITH COV 0.5
    pp.check.hab.cov0.5 <- 
      ggarrange(get_ppcheck("Complex ~ NePRA binary", m.c.neprab.hab.cov0.5),
                get_ppcheck("Simple ~ NePRA binary", m.c.neprab.hab.cov0.5),
                get_ppcheck("Non-tool ~ NePRA binary", m.c.neprab.hab.cov0.5),
                get_ppcheck("Complex ~ NePRA similarity", m.c.neprasim.hab.cov0.5),
                get_ppcheck("Simple ~ NePRA similarity", m.c.neprasim.hab.cov0.5),
                get_ppcheck("Non-tool ~ NePRA similarity", m.c.neprasim.hab.cov0.5),
                get_ppcheck("Complex ~ IBD binary", m.c.ibd.hab.cov0.5),
                get_ppcheck("Simple ~ IBD binary", m.c.ibd.hab.cov0.5),
                get_ppcheck("Non-tool ~ IBD binary", m.c.ibd.hab.cov0.5),
                nrow = 3, ncol = 3)
    
    ggsave(plot = pp.check.hab.cov0.5, file = "Results/Bayesian/md_ppcheck_habitat_cov0.5.jpeg", 
           dpi = 300, width = 15, height = 15, unit = "cm")   

  # MODELS WITH COV 0.1
    pp.check.hab.cov0.1 <- 
      ggarrange(get_ppcheck("Complex ~ NePRA binary", m.c.neprab.hab.cov0.1),
                get_ppcheck("Simple ~ NePRA binary", m.c.neprab.hab.cov0.1),
                get_ppcheck("Non-tool ~ NePRA binary", m.c.neprab.hab.cov0.1),
                get_ppcheck("Complex ~ NePRA similarity", m.c.neprasim.hab.cov0.1),
                get_ppcheck("Simple ~ NePRA similarity", m.c.neprasim.hab.cov0.1),
                get_ppcheck("Non-tool ~ NePRA similarity", m.c.neprasim.hab.cov0.1),
                get_ppcheck("Complex ~ IBD binary", m.c.ibd.hab.cov0.1),
                get_ppcheck("Simple ~ IBD binary", m.c.ibd.hab.cov0.1),
                get_ppcheck("Non-tool ~ IBD binary", m.c.ibd.hab.cov0.1),
                nrow = 3, ncol = 3)
    
    ggsave(plot = pp.check.hab.cov0.5, file = "Results/Bayesian/md_ppcheck_habitat_cov0.1.jpeg", 
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


  # GET PLOTS AND SAVE AS PDF OF COV 0.5
    pdf("Results/Bayesian/md_tracedensity_habitat_cov0.5.pdf", width = 15, height = 10)
    
    
    # Arrange pages
    get_plots("Complex ~ NePRA binary", m.c.neprab.hab.cov0.5,
              c("b_Intercept", "b_neprabinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Simple ~ NePRA binary", m.s.neprab.hab.cov0.5,
              c("b_Intercept", "b_neprabinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Non-tool ~ NePRA binary", m.n.neprab.hab.cov0.5,
              c("b_Intercept", "b_neprabinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Complex ~ NePRA similarity", m.c.neprasim.hab.cov0.5,
              c("b_Intercept", "b_z.log.nepra", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Simple ~ NePRA similarity", m.s.neprasim.hab.cov0.5,
              c("b_Intercept", "b_z.log.nepra", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept")) 
    get_plots("Non-tool ~ NePRA similarity", m.n.neprasim.hab.cov0.5,
              c("b_Intercept", "b_z.log.nepra", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Complex ~ IBD binary", m.c.ibd.hab.cov0.5,
              c("b_Intercept", "b_ibdbinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Simple ~ IBD binary", m.s.ibd.hab.cov0.5,
              c("b_Intercept", "b_ibdbinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept")) 
    get_plots("Non-tool ~ IBD binary", m.n.ibd.hab.cov0.5,
              c("b_Intercept", "b_ibdbinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept")) 
    dev.off()
    
  # GET PLOTS AND SAVE AS PDF OF COV 0.1
    pdf("Results/Bayesian/md_tracedensity_habitat_cov0.1.pdf", width = 15, height = 10)
    
    
    # Arrange pages
    get_plots("Complex ~ NePRA binary", m.c.neprab.hab.cov0.1,
              c("b_Intercept", "b_neprabinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Simple ~ NePRA binary", m.s.neprab.hab.cov0.1,
              c("b_Intercept", "b_neprabinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Non-tool ~ NePRA binary", m.n.neprab.hab.cov0.1,
              c("b_Intercept", "b_neprabinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Complex ~ NePRA similarity", m.c.neprasim.hab.cov0.1,
              c("b_Intercept", "b_z.log.nepra", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Simple ~ NePRA similarity", m.s.neprasim.hab.cov0.1,
              c("b_Intercept", "b_z.log.nepra", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept")) 
    get_plots("Non-tool ~ NePRA similarity", m.n.neprasim.hab.cov0.1,
              c("b_Intercept", "b_z.log.nepra", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Complex ~ IBD binary", m.c.ibd.hab.cov0.1,
              c("b_Intercept", "b_ibdbinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept"))
    get_plots("Simple ~ IBD binary", m.s.ibd.hab.cov0.1,
              c("b_Intercept", "b_ibdbinary", "b_samehabitat", "b_z.obstime", "sd_dyadID__Intercept", "sd_behaviour__Intercept")) 
    get_plots("Non-tool ~ IBD binary", m.n.ibd.hab.cov0.1,
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
  
  # GET PLOTS WITH COV 0.5

    pdf("Results/Bayesian/md_loo_habitat_cov0.5.pdf", width = 8, height = 7)
    
    par(mfrow=c(3,3))
    
    get_pareto("Complex ~ NePRA binary", m.c.neprab.hab.cov0.5)
    get_pareto("Complex ~ NePRA similarity", m.c.neprasim.hab.cov0.5)
    get_pareto("Complex ~ IBD binary", m.c.ibd.hab.cov0.5)
    get_pareto("Simple ~ NePRA binary", m.s.neprab.hab.cov0.5)
    get_pareto("Simple ~ NePRA similarity", m.s.neprasim.hab.cov0.5)
    get_pareto("Simple ~ IBD binary", m.s.ibd.hab.cov0.5)
    get_pareto("Non-tool ~ NePRA binary", m.n.neprab.hab.cov0.5)
    get_pareto("Non-tool ~ NePRA similarity", m.n.neprasim.hab.cov0.5)
    get_pareto("Non-tool ~ IBD binary", m.n.ibd.hab.cov0.5)
    
    dev.off()
    
    
  # GET PLOTS WITH COV 0.1
    
    pdf("Results/Bayesian/md_loo_habitat_cov0.1.pdf", width = 8, height = 7)
    
    par(mfrow=c(3,3))
    
    get_pareto("Complex ~ NePRA binary", m.c.neprab.hab.cov0.1)
    get_pareto("Complex ~ NePRA similarity", m.c.neprasim.hab.cov0.1)
    get_pareto("Complex ~ IBD binary", m.c.ibd.hab.cov0.1)
    get_pareto("Simple ~ NePRA binary", m.s.neprab.hab.cov0.1)
    get_pareto("Simple ~ NePRA similarity", m.s.neprasim.hab.cov0.1)
    get_pareto("Simple ~ IBD binary", m.s.ibd.hab.cov0.1)
    get_pareto("Non-tool ~ NePRA binary", m.n.neprab.hab.cov0.1)
    get_pareto("Non-tool ~ NePRA similarity", m.n.neprasim.hab.cov0.1)
    get_pareto("Non-tool ~ IBD binary", m.n.ibd.hab.cov0.1)
    
    dev.off()
    
############################################### MODEL COMPARISON WITH LOO & WAIC -----
    
  # IMPORT ORIGINAL MODELS
    load("Results/Bayesian/models_habitat.RData")

  # LOO AND WAIC VALUES
    loo(m.c.neprab.hab)
    waic(m.c.neprab.hab)
    
    loo(m.c.neprab.hab.cov0.5)
    waic(m.c.neprab.hab.cov0.5)
    
    loo(m.c.neprab.hab.cov0.1)
    waic(m.c.neprab.hab.cov0.1)
    
    
    loo(m.c.neprasim.hab)
    waic(m.c.neprasim.hab)
    
    loo(m.c.neprasim.hab.cov0.5)
    waic(m.c.neprasim.hab.cov0.5)
    
    loo(m.c.neprasim.hab.cov0.1)
    waic(m.c.neprasim.hab.cov0.1)
    
    loo(m.c.ibd.hab)
    waic(m.c.ibd.hab)
    
    loo(m.c.ibd.hab.cov0.5)
    waic(m.c.ibd.hab.cov0.5)
    
    loo(m.c.ibd.hab.cov0.1)
    waic(m.c.ibd.hab.cov0.1)
    
  # DIRECT COMPARISON
    m.c.neprab.hab <- add_criterion(m.c.neprab.hab, "loo")
    m.c.neprab.hab.cov0.5 <- add_criterion(m.c.neprab.hab.cov0.5, "loo")
    m.c.neprab.hab.cov0.1 <- add_criterion(m.c.neprab.hab.cov0.1, "loo")

    loo_compare(m.c.neprab.hab, m.c.neprab.hab.cov0.5, m.c.neprab.hab.cov0.1,
                criterion = "loo")
    

