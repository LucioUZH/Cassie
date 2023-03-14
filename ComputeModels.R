# ==============================================================================
# Cassandra Gunasekaram, AIM UZH, PhD - Chimp genetics project, 21.02.2023
# Run tests/models
# Last update: 14.03.2023
# ==============================================================================

# ______________________________________________________________________________
# Workspace preparations and data input----

  rm(list=ls())
  setwd("C:/Users/cassi/OneDrive - Universität Zürich UZH/PhD/Genetics/Analysis/Repository")
  library(dplyr)
  library(ggpubr)
  library(cgwtools) # too resave

## GENERAL FUNCTIONS............................................................
  # Not in function
  '%!in%' <- function(x,y)!('%in%'(x,y))

## IMPORT DATA..................................................................
  ddyads <- select(read.csv("Data/d.dyads.csv"), -1)
  ddyads.raw <- ddyads

# ----
# ______________________________________________________________________________
# Bayesian models-----

## DATA PREPARATIONS............................................................

# In ecology paper, they exclude non-applicable cases - I did not do that - consider this?
  
# Tidy and create parameters
  ddyads$behaviour <- as.factor(ddyads$behaviour)
  ddyads$category <- as.factor(ddyads$category)
  ddyads$z.totalobstime <- scale(ddyads$totalobstime)
  ddyads$z.log.similarity <- scale(log(ddyads$similarity))
  ddyads$z.sqrt.meanibd <- scale(sqrt(ddyads$meanibd)) # change this to z.log
  ddyads$z.sqrt.maxibdc <- scale(sqrt(ddyads$maxIBDc)) # change this to avg.maxibd and z.log
  
  ddyads$sameregion <- factor(ddyads$sameregion, levels = c("1", "0"))
  ddyads$samehabitat <- factor(ddyads$samehabitat, levels = c("1", "0"))

# Subset data to nontool, tools, simple and complex tools
  dntool <- droplevels(filter(ddyads, complexity == 0))
  dtools <- droplevels(filter(ddyads, complexity > 0))
  dsimple <- droplevels(filter(ddyads, complexity == 1))
  dcomplex <- droplevels(filter(ddyads, complexity == 2))

## BAYESIAN MODEL PREPARATIONS..................................................

  library(brms)
  library(rstan)
  library(tidybayes)

# In case it crashes
  rstan_options(auto_write = TRUE)

# Function to set prior
  setprior <- function(dist,testvar){
    priors <- c(set_prior(dist, class = "Intercept"),
               set_prior(dist, class = "b"),
               set_prior(dist, class = "b", coef = testvar),
               set_prior(dist, class = "b", coef = "z.totalobstime"),
               set_prior(dist, class="sd")
    )
    return(priors)
  }
  # hd.prior <- c(prior(normal(0,1), class=Intercept),
  #               prior(normal(0,1), class=b),
  #               prior(exponential(1), class=sd),
  #               prior(normal(0,1), class=Intercept, dpar="hu"),
  #               prior(normal(0,1), class=b, dpar="hu"),
  #               prior(exponential(1), class=sd, dpar="hu"))


# Function to run bayesian model
  bayesianmodel <- function(formel, fam, dataset, priors){
    bm <- brm(formula = formel,family = fam,data = dataset,
              prior= priors, 
              cores = 16, seed = 100, 
              warmup=1000, iter=2000, chains=4,
              control=list(adapt_delta =0.99))
  }

## MODELS TESTING RARE ALLELE SIMILARITY.......................................

# Function
  f <- bf(shared ~ z.log.similarity + 
            z.totalobstime + sameregion + samehabitat + (1|dyadID) + 
            (1+ z.log.similarity + z.totalobstime + 
               dyadregion + samehabitat|behaviour))

# Default prior
  priors <- setprior("normal(0,1)", "z.log.similarity")

  bm.ntool.snp <- bayesianmodel(f, "bernoulli",dntool,priors)
  bm.tool.snp <- bayesianmodel(f, "bernoulli",dtool,priors)
  bm.simple.snp <- bayesianmodel(f, "bernoulli",dsimple,priors)
  bm.complex.snp <- bayesianmodel(f, "bernoulli",dcomplex,priors)

  # save(bm.ntool.snp, bm.tool.snp, bm.simple.snp, bm.complex.snp,
  #      file = "Results/bayesiandefault.RData")
  remove(priors, bm.ntool.snp,bm.tool.snp,bm.simple.snp,bm.complex.snp)

# Wide prior
  priors <- setprior("normal(0,10)", "z.log.similarity")
  
  bmwide.ntool.snp <- bayesianmodel(f, "bernoulli",dntool,priors)
  bmwide.tool.snp <- bayesianmodel(f, "bernoulli",dtools,priors)
  bmwide.simple.snp <- bayesianmodel(f, "bernoulli",dsimple,priors)
  bmwide.complex.snp <- bayesianmodel(f, "bernoulli",dcomplex,priors)

  # save(bmwide.ntool.snp, bmwide.tool.snp, bmwide.simple.snp, bmwide.complex.snp,
  #      file = "Results/bayesianwide.RData")
  remove(f, priors, bmwide.ntool.snp, bmwide.tool.snp, bmwide.simple.snp, bmwide.complex.snp)


## MODELS TESTING IBD Y/N.......................................................

# Function
  f <- bf(shared ~ presenceibd + 
            z.totalobstime + dyadregion + samehabitat + (1|dyadID) + 
            (1+ presenceibd + z.totalobstime + 
               dyadregion + samehabitat|behaviour))

# Default prior
  priors <- setprior("normal(0,1)", "presenceibd")

  bm.ntool.ibdyn <- bayesianmodel(f, "bernoulli",dntool,priors)
  bm.tool.ibdyn <- bayesianmodel(f, "bernoulli",dtools,priors)
  bm.simple.ibdyn <- bayesianmodel(f, "bernoulli",dsimple,priors)
  bm.complex.ibdyn <- bayesianmodel(f, "bernoulli",dcomplex,priors)

  resave(bm.ntool.ibdyn, bm.tool.ibdyn, bm.simple.ibdyn, bm.complex.ibdyn,
         file = "Results/bayesiandefault.RData")
  remove(priors,bm.ntool.ibdyn,bm.tool.ibdyn,bm.simple.ibdyn,bm.complex.ibdyn)

# Wide prior
  priors <- setprior("normal(0,10)", "presenceibd")
  
  bmwide.ntool.ibdyn <- bayesianmodel(f, "bernoulli",dntool,priors)
  bmwide.tool.ibdyn <- bayesianmodel(f, "bernoulli",dtools,priors)
  bmwide.simple.ibdyn <- bayesianmodel(f, "bernoulli",dsimple,priors)
  bmwide.complex.ibdyn <- bayesianmodel(f, "bernoulli",dcomplex,priors)
  
  resave(bmwide.ntool.ibdyn,bmwide.tool.ibdyn,bmwide.simple.ibdyn,bmwide.complex.ibdyn,
         file = "Results/bayesianwide.RData")
  remove(f,priors,bmwide.ntool.ibdyn,bmwide.tool.ibdyn,bmwide.simple.ibdyn,bmwide.complex.ibdyn)

## MODELS TESTING MEAN IBD Length.......................................................

# Function
  f <- bf(shared ~ z.sqrt.meanibd + 
            z.totalobstime + samehabitat + (1|dyadID) + 
            (z.sqrt.meanibd|dyadregion)+
            (1+ z.sqrt.meanibd + z.totalobstime + 
               samehabitat|behaviour))

# Default prior
  priors <- setprior("normal(0,1)", "z.sqrt.meanibd")

  #! subset data to presenceibd==1
  bm.ntool.ibdmean <- bayesianmodel(f, "bernoulli",dntool,priors)
  bm.tool.ibdmean <- bayesianmodel(f, "bernoulli",dtools,priors)
  bm.simple.ibdmean <- bayesianmodel(f, "bernoulli",dsimple,priors)
  bm.complex.ibdmean <- bayesianmodel(f, "bernoulli",dcomplex,priors)

  resave(bm.ntool.ibdmean, bm.tool.ibdmean, bm.simple.ibdmean, bm.complex.ibdmean,
         file = "Results/bayesiandefault.RData")
  remove(priors,bm.ntool.ibdmean, bm.tool.ibdmean, bm.simple.ibdmean, bm.complex.ibdmean)

# Wide prior
  priors <- setprior("normal(0,10)", "z.sqrt.meanibd")
  
  #! again subset data to presenceibd==1
  bmwide.ntool.ibdmean <- bayesianmodel(f, "bernoulli",dntool,priors)
  bmwide.tool.ibdmean <- bayesianmodel(f, "bernoulli",dtools,priors)
  bmwide.simple.ibdmean <- bayesianmodel(f, "bernoulli",dsimple,priors)
  bmwide.complex.ibdmean <- bayesianmodel(f, "bernoulli",dcomplex,priors)

  resave(bmwide.ntool.ibdmean, bmwide.tool.ibdmean, bmwide.simple.ibdmean, bmwide.complex.ibdmean,
         file = "Results/bayesianwide.RData")
  remove(f,priors,bmwide.ntool.ibdmean, bmwide.tool.ibdmean, bmwide.simple.ibdmean, bmwide.complex.ibdmean)

## MODELS TESTING Max IBD Length.......................................................

  #! change this to avg.maxibd
# Function
  f <- bf(shared ~ z.sqrt.maxibdc + 
            z.totalobstime + dyadregion + samehabitat + (1|dyadID) + 
            (1+ z.sqrt.maxibdc + z.totalobstime + 
               dyadregion + samehabitat|behaviour))

# Default prior
  priors <- setprior("normal(0,1)", "z.sqrt.maxibd")
  
  bm.ntool.ibdmax <- bayesianmodel(f, "bernoulli",dntool,priors)
  bm.tool.ibdmax <- bayesianmodel(f, "bernoulli",dtools,priors)
  bm.simple.ibdmax <- bayesianmodel(f, "bernoulli",dsimple,priors)
  bm.complex.ibdmax <- bayesianmodel(f, "bernoulli",dcomplex,priors)

  resave(bm.ntool.ibdmax, bm.tool.ibdmax, bm.simple.ibdmax, bm.complex.ibdmax,
         file = "Results/bayesiandefault.RData")
  remove(priors,bm.ntool.ibdmax, bm.tool.ibdmax, bm.simple.ibdmax, bm.complex.ibdmax)

# Wide prior
  priors <- setprior("normal(0,10)", "z.sqrt.meanibd")
  
  bmwide.ntool.ibdmax <- bayesianmodel(f, "bernoulli",dntool,priors)
  bmwide.tool.ibdmax <- bayesianmodel(f, "bernoulli",dtools,priors)
  bmwide.simple.ibdmax <- bayesianmodel(f, "bernoulli",dsimple,priors)
  bmwide.complex.ibdmax <- bayesianmodel(f, "bernoulli",dcomplex,priors)
  
  resave(bmwide.ntool.ibdmax, bmwide.tool.ibdmax, bmwide.simple.ibdmax, bmwide.complex.ibdmax,
         file = "Results/bayesianwide.RData")
  remove(f,priors,bmwide.ntool.ibdmax, bmwide.tool.ibdmax, bmwide.simple.ibdmax, bmwide.complex.ibdmax)

## CLEAR WORKSPACE..............................................................
  rm(list= ls()[!(ls() %in% c('ddyads.raw','ddyads', "%!in%"))])

  detach("package:brms", unload=TRUE)
  detach("package:rstan", unload=TRUE)
  detach("package:tidybayes", unload=TRUE)
  


# ----
# ______________________________________________________________________________
# Mantel correlation tests & correlogram----

## PREPARE MATRICES.............................................................
# Function to create genetic and geographic matrices
  creatematrix <- function(dd, variable, dpop){
    dd <- distinct(dd, dyadID,.keep_all = T)
    dd <- merge(dd, dpop, by.x = "population1", by.y = "pop")
    dd <- merge(dd, dpop, by.x = "population2", by.y = "pop")
    dd2 <- select(dd, c("population1", "population2", "index.x", 
                        "index.y", "dyadID", variable))
    colnames(dd2)[1:4] <- c("population2", "population1", "index.y", "index.x")
    dd2 <- rbind(select(dd, c("population1", "population2", "index.x", 
                              "index.y", "dyadID", variable)), dd2)
    output <- data.frame(matrix(NA, nrow = as.numeric(nrow(dpop)), 
                                ncol = as.numeric(nrow(dpop))))
  
    for (i in 1:nrow(dd2)) {
      x <- dd2$index.x[i]
      y <- dd2$index.y[i]
      output[dd2$index.x[i],dd2$index.y[i]] <- dd2[i,6]
    }
    colnames(output) <- dpop$pop
    return(output)
  }

# Prepare indices for populations to match with each other
  pops <- sort(unique(c(ddyads$population1, ddyads$population2)))
  pops <- data.frame(pop = pops, index = seq(1,34))

# Nontool dataframes
  mibd <- creatematrix(ddyads, "avgmaxibd",pops)
  msimilarity <- creatematrix(ddyads, "similarity", pops)
  mgeo <- creatematrix(ddyads, "distance", pops)
  mregion <- creatematrix(ddyads, "sameregion", pops)
  mhab <- creatematrix(ddyads, "samehabitat", pops)

# Import cultural data
  dcult.raw <- read.csv("Data/d.culture.csv", sep = ";")
  
# Function to create complexity cultural matrices
  creatematrix2 <- function(data, clevel){
    x <- droplevels(filter(data, complexity == clevel))
    nbehav <- nrow(distinct(x, behaviour))
    mcult <- data.frame(matrix(NA, nrow = 34, ncol = nbehav))
    colnames(mcult) <- distinct(x,behaviour)$behaviour
    mcult$pop <- distinct(x,site)$site
    for (i in 1:nrow(x)){
      toolindex <- which(colnames(mcult) == x$behaviour[i])
      popindex <- which(mcult$pop == x$site[i])
      mcult[popindex,toolindex] <- x$presencewithout[i]
    }
    mcult <- select(mcult, -last_col())
    rownames(mcult) <- distinct(x, site)$site
    return(mcult)
    remove(x,nbehav,mcult,i,toolindex,popindex)
  }

# Create complexity cultural data
  mcult0 <- creatematrix2(dcult.raw,0)
  mcult1 <- creatematrix2(dcult.raw,1)
  mcult2 <- creatematrix2(dcult.raw,2)

# Function to create complex  matrices
  creatematrix3 <- function(data, clevel, cat){
    x <- droplevels(filter(data, complexity == clevel & category %in% cat))
    nbehav <- nrow(distinct(x, behaviour))
    mcult <- data.frame(matrix(NA, nrow = 34, ncol = nbehav))
    colnames(mcult) <- distinct(x,behaviour)$behaviour
    mcult$pop <- distinct(x,site)$site
    for (i in 1:nrow(x)){
      toolindex <- which(colnames(mcult) == x$behaviour[i])
      popindex <- which(mcult$pop == x$site[i])
      mcult[popindex,toolindex] <- x$presencewithout[i]
    }
    mcult <- select(mcult, -last_col())
    rownames(mcult) <- distinct(x, site)$site
    return(mcult)
    remove(x,nbehav,mcult,i,toolindex,popindex)
  }

# Create complex dataframes
  mcultu <- creatematrix3(dcult.raw,2, c("perforate honey/termites"))
  mcultn <- creatematrix3(dcult.raw,2, c("nut"))
  mculth <- creatematrix3(dcult.raw,2, c("honey"))
  mcultc <- creatematrix3(dcult.raw,2, c("communication"))

# Simple simple dataframes
  mcultsc <- creatematrix3(dcult.raw,1, c("communication"))
  mcultsi <- creatematrix3(dcult.raw,1, c("ant", "honey", "termite"))
  mcultsf <- creatematrix3(dcult.raw,1, c("algae", "water", "vertebrates"))


## DISTANCE MATRICES...................................................

  library(vegan)

# Genetic and geographic matrices
  simdist <- as.dist(msimilarity, diag=F, upper=F)
  ibddist <- as.dist(mibd, diag=F, upper=F)
  geodist <- as.dist(mgeo, diag=F, upper=F)
  regiondist <- as.dist(mregion, diag=F, upper=F)
  habdist <- as.dist(mhab, diag=F, upper=F)

# Function for cultural matrices implementing Jaccard disimilarity matrix
  cultdist <- function(m){
    jacdist <- vegdist(m, method = "jaccard", diag=F, upper=F)
    # set all NaN (0/0) to 1.0 (maximum dissimilarity)
    jacdist[is.nan(jacdist)] = 1.0
    # for more intuitive interpretation of results, calculate similarity from dissimilarity
    jacdist <- 1-jacdist
  }

# Jaccard dissimliarity matrix for cultural traits
  cultdist0 <- cultdist(mcult0)
  cultdist1 <- cultdist(mcult1)
  cultdist2 <- cultdist(mcult2)
  
  cultdistc <- cultdist(mcultc)
  cultdisth <- cultdist(mculth)
  cultdistn <- cultdist(mcultn)
  cultdistu <- cultdist(mcultu)
  
  cultdistsc <- cultdist(mcultsc)
  cultdistsi <- cultdist(mcultsi)
  cultdistsf <- cultdist(mcultsf)

# Save distances in data frame 
  distframe <- data.frame('Euclidian\ndistance'=as.vector(geodist),
                          'Same region'=as.vector(regiondist),
                          'Same habitat'=as.vector(habdist),
                          'Rare alleles' =as.vector(simdist),
                          'Mean IBD\nlength'=as.vector(ibddist),
                          'Non-tool\nbehaviours'=as.vector(cultdist0),
                          'Simple tools'=as.vector(cultdist1),
                          'Composite tools'=as.vector(cultdist2),
                          'Acc. stone\nthrowing'=as.vector(cultdistc),
                          'Honey toolset'=as.vector(cultdisth),
                          'Nut cracking'=as.vector(cultdistn),
                          'Underground'=as.vector(cultdistu)
  )

  save(geodist, habdist, regiondist, ibddist, simdist, cultdist0, cultdist1, cultdist2,
       cultdisth, cultdistn, cultdistc, cultdistu, cultdistsf, cultdistsi, cultdistsc,
       file = "Results/similaritymatrices.RData")
## SIMPLE MANTEL TEST...........................................................

  detach("package:vegan", unload=TRUE)
  library(ecodist)

# Set up matrices to store results # This code bit comes from the pottery paper
  rows = 7
  cols = 7
  cor_colnames<-c('Euclidean\ndistance',
                  'Habitat\n(Same/different)','Rare alleles', 
                  'Avg. Maximum\nlength','Non-tool\nbehaviours',
                  'Simple tools','Composite tools')
  
  mantel.all.r<-matrix(1,nrow=rows,ncol=cols)
  rownames(mantel.all.r)<-cor_colnames
  colnames(mantel.all.r)<-cor_colnames
  mantel.all.cil<-mantel.all.ciu<-mantel.all.p<-mantel.all.r
  
  store_result<-function(i,j,MTR) {
    mantel.all.cil[i,j] <<- MTR[5]
    mantel.all.ciu[i,j] <<- MTR[6]
    mantel.all.p[i,j] <<- MTR[4]
    mantel.all.r[i,j] <<- MTR[1]
  }

# Function to run through correlation tests and store data in the prepared matrices
  mantelmodel <- function(distmatrices){
    for (i in 1:6){
      n=i
      dist1 <- distmatrices[[i]]
      while (n<7){
        n = n+1
        dist2 <- distmatrices[[n]]
        test <- mantel(dist1 ~ dist2, nperm = 10000, mrank = T)
        store_result(n,i,test)
      }
    }
  }

# Run Mantel tests
  matrices <- list(geodist, habdist, simdist, ibddist, cultdist0, cultdist1, cultdist2)
  mantelmodel(matrices)
  for(i in 1:rows) for(j in 1:cols) mantel.all.r[i,j]<-mantel.all.r[j,i]
  for(i in 1:rows) for(j in 1:cols) mantel.all.p[i,j]<-mantel.all.p[j,i]
  
# Save and remove
  # save(mantel.all.cil, mantel.all.ciu, mantel.all.p, mantel.all.r,
  #      file = "Results/manteltests.RData")
  remove(mantel.all.cil,mantel.all.ciu,mantel.all.p,mantel.all.r)

## PARTIAL MANTEL TESTS.........................................................

# IBD
  rows = 5
  cols = 5
  Pcor_colnames<-c('Region\n(same/different)', "Rare alleles",
                   'Non-tool\nbehaviours','Simple tools',
                   'Composite tools')

  mantel.ibd.r<-matrix(1,nrow=rows,ncol=cols)
  rownames(mantel.ibd.r)<-Pcor_colnames
  colnames(mantel.ibd.r)<-Pcor_colnames
  mantel.ibd.cil<-mantel.ibd.ciu<-mantel.ibd.p<-mantel.ibd.r

  store_result2<-function(i,j,MTR) {
    mantel.ibd.cil[i,j] <<- MTR[5]
    mantel.ibd.ciu[i,j] <<- MTR[6]
    mantel.ibd.p[i,j] <<- MTR[4]
    mantel.ibd.r[i,j] <<- MTR[1]
  }

# Function to run through correlation tests and store data in the prepared matrices
  ibdmantelmodel <- function(distmatrices, genedist){
    for (i in 1:4){
      n=i
      dist1 <- distmatrices[[i]]
      while (n<5){
        n = n+1
        dist2 <- distmatrices[[n]]
        test <- mantel(dist1 ~ dist2 + genedist, nperm = 10000, mrank = T)
        store_result2(n,i,test)
      }
    }
  }

# Run partial Mantel test
  matrices <- list(regiondist, simdist, cultdist0, cultdist1, cultdist2)
  ibdmantelmodel(matrices, ibddist)
  for(i in 1:rows) for(j in 1:cols) mantel.ibd.r[i,j]<-mantel.ibd.r[j,i]
  for(i in 1:rows) for(j in 1:cols) mantel.ibd.p[i,j]<-mantel.ibd.p[j,i]

# Save and remove
  resave(mantel.ibd.cil, mantel.ibd.ciu, mantel.ibd.p, mantel.ibd.r,
       file = "Results/manteltests.RData")
  remove(mantel.ibd.cil, mantel.ibd.ciu, mantel.ibd.p, mantel.ibd.r)

# Same for rare alleles
  Pcor_colnames<-c('Region\n(same/different)','Mean IBD length',
                   'Non-tool\nbehaviours','Simple tools',
                   'Composite tools')
  mantel.snp.r<-matrix(1,nrow=rows,ncol=cols)
  rownames(mantel.snp.r)<-Pcor_colnames
  colnames(mantel.snp.r)<-Pcor_colnames
  mantel.snp.cil<-mantel.snp.ciu<-mantel.snp.p<-mantel.snp.r
  
  store_result3<-function(i,j,MTR) {
    mantel.snp.cil[i,j] <<- MTR[5]
    mantel.snp.ciu[i,j] <<- MTR[6]
    mantel.snp.p[i,j] <<- MTR[4]
    mantel.snp.r[i,j] <<- MTR[1]
  }
  
  snpmantelmodel <- function(distmatrices, genedist){
    for (i in 1:4){
      n=i
      dist1 <- distmatrices[[i]]
      while (n<5){
        n = n+1
        dist2 <- distmatrices[[n]]
        test <- mantel(dist1 ~ dist2 + genedist, nperm = 10000, mrank = T)
        store_result3(n,i,test)
      }
    }
  }

# Run partial Mantel test
  matrices <- list(regiondist, ibddist, cultdist0, cultdist1, cultdist2)
  snpmantelmodel(matrices, simdist)
  for(i in 1:rows) for(j in 1:cols) mantel.snp.r[i,j]<-mantel.snp.r[j,i]
  for(i in 1:rows) for(j in 1:cols) mantel.snp.p[i,j]<-mantel.snp.p[j,i]
  
# Save and remove
  resave(mantel.snp.cil, mantel.snp.ciu, mantel.snp.p, mantel.snp.r,
         file = "Results/manteltests.RData")
  remove(mantel.snp.cil, mantel.snp.ciu, mantel.snp.p, mantel.snp.r)
  

## SIMPLE MANTEL WITH COMPLEX TOOLS.............................................
  rows = 3
  cols = 3

# Acc. stone throwing
  cor_colnames<-c('Rare alleles', 
                  'Avg. maximum\nIBD length','Acc. stone\nthrowing')
  mantel.accstone.r<-matrix(1,nrow=rows,ncol=cols)
  rownames(mantel.accstone.r)<-cor_colnames
  colnames(mantel.accstone.r)<-cor_colnames
  mantel.accstone.cil<-mantel.accstone.ciu<-mantel.accstone.p<-mantel.accstone.r

  store_resultC<-function(i,j,MTR) {
    mantel.accstone.cil[i,j] <<- MTR[5]
    mantel.accstone.ciu[i,j] <<- MTR[6]
    mantel.accstone.p[i,j] <<- MTR[4]
    mantel.accstone.r[i,j] <<- MTR[1]
  }

  store_resultC(2,1,mantel(simdist ~ ibddist, nperm = 10000, mrank = T))
  store_resultC(3,1,mantel(simdist ~ cultdistc, nperm = 10000, mrank = T))
  store_resultC(3,2,mantel(ibddist ~ cultdistc, nperm = 10000, mrank = T))

  for(i in 1:rows) for(j in 1:cols) mantel.accstone.r[i,j]<-mantel.accstone.r[j,i]
  for(i in 1:rows) for(j in 1:cols) mantel.accstone.p[i,j]<-mantel.accstone.p[j,i]
  
# Save and remove
  resave(mantel.accstone.cil, mantel.accstone.ciu, mantel.accstone.p, mantel.accstone.r,
         file = "Results/manteltests.RData")
  remove(mantel.accstone.cil, mantel.accstone.ciu, mantel.accstone.p, mantel.accstone.r)
  

# Honey toolset
  cor_colnames<-c('Rare alleles', 
                  'Avg. maximum\nIBD length','Honey toolset')
  mantel.honey.r<-matrix(1,nrow=rows,ncol=cols)
  rownames(mantel.honey.r)<-cor_colnames
  colnames(mantel.honey.r)<-cor_colnames
  mantel.honey.cil<-mantel.honey.ciu<-mantel.honey.p<-mantel.honey.r
  
  store_resultC<-function(i,j,MTR) {
    mantel.honey.cil[i,j] <<- MTR[5]
    mantel.honey.ciu[i,j] <<- MTR[6]
    mantel.honey.p[i,j] <<- MTR[4]
    mantel.honey.r[i,j] <<- MTR[1]
  }

  store_resultC(2,1,mantel(simdist ~ ibddist, nperm = 10000, mrank = T))
  store_resultC(3,1,mantel(simdist ~ cultdisth, nperm = 10000, mrank = T))
  store_resultC(3,2,mantel(ibddist ~ cultdisth, nperm = 10000, mrank = T))

  for(i in 1:rows) for(j in 1:cols) mantel.honey.r[i,j]<-mantel.honey.r[j,i]
  for(i in 1:rows) for(j in 1:cols) mantel.honey.p[i,j]<-mantel.honey.p[j,i]

# Save and remove
  resave(mantel.honey.cil, mantel.honey.ciu, mantel.honey.p, mantel.honey.r,
         file = "Results/manteltests.RData")
  remove(mantel.honey.cil, mantel.honey.ciu, mantel.honey.p, mantel.honey.r)
  

# Nut cracking 
  cor_colnames<-c('Rare alleles', 
                  'Avg. maximum\nIBD length','Nut cracking')
  mantel.nut.r<-matrix(1,nrow=rows,ncol=cols)
  rownames(mantel.nut.r)<-cor_colnames
  colnames(mantel.nut.r)<-cor_colnames
  mantel.nut.cil<-mantel.nut.ciu<-mantel.nut.p<-mantel.nut.r

  store_resultC<-function(i,j,MTR) {
    mantel.nut.cil[i,j] <<- MTR[5]
    mantel.nut.ciu[i,j] <<- MTR[6]
    mantel.nut.p[i,j] <<- MTR[4]
    mantel.nut.r[i,j] <<- MTR[1]
  }

  store_resultC(2,1,mantel(simdist ~ ibddist, nperm = 10000, mrank = T))
  store_resultC(3,1,mantel(simdist ~ cultdistn, nperm = 10000, mrank = T))
  store_resultC(3,2,mantel(ibddist ~ cultdistn, nperm = 10000, mrank = T))

  for(i in 1:rows) for(j in 1:cols) mantel.nut.r[i,j]<-mantel.nut.r[j,i]
  for(i in 1:rows) for(j in 1:cols) mantel.nut.p[i,j]<-mantel.nut.p[j,i]
  
# Save and remove
  resave(mantel.nut.cil, mantel.nut.ciu, mantel.nut.p, mantel.nut.r,
         file = "Results/manteltests.RData")
  remove(mantel.nut.cil, mantel.nut.ciu, mantel.nut.p, mantel.nut.r)
  

# Underground
  cor_colnames<-c('Rare alleles', 
                  'Avg. maximum\nIBD length','Underground')
  mantel.underground.r<-matrix(1,nrow=rows,ncol=cols)
  rownames(mantel.underground.r)<-cor_colnames
  colnames(mantel.underground.r)<-cor_colnames
  mantel.underground.cil<-mantel.underground.ciu<-mantel.underground.p<-mantel.underground.r

  store_resultC<-function(i,j,MTR) {
    mantel.underground.cil[i,j] <<- MTR[5]
    mantel.underground.ciu[i,j] <<- MTR[6]
    mantel.underground.p[i,j] <<- MTR[4]
    mantel.underground.r[i,j] <<- MTR[1]
  }

  store_resultC(2,1,mantel(simdist ~ ibddist, nperm = 10000, mrank = T))
  store_resultC(3,1,mantel(simdist ~ cultdistu, nperm = 10000, mrank = T))
  store_resultC(3,2,mantel(ibddist ~ cultdistu, nperm = 10000, mrank = T))

  for(i in 1:rows) for(j in 1:cols) mantel.underground.r[i,j]<-mantel.underground.r[j,i]
  for(i in 1:rows) for(j in 1:cols) mantel.underground.p[i,j]<-mantel.underground.p[j,i]
  
# Save and remove
  resave(mantel.underground.cil, mantel.underground.ciu, mantel.underground.p, mantel.underground.r,
         file = "Results/manteltests.RData")
  remove(mantel.underground.cil, mantel.underground.ciu, mantel.underground.p, mantel.underground.r)
  

## SIMPLE MANTEL WITH SIMPLE TOOLS.............................................
  rows = 3
  cols = 3
  
# Communication
  cor_colnames<-c('Rare alleles', 
                  'Mean IBD length','Communication')
  mantel.comm.r<-matrix(1,nrow=rows,ncol=cols)
  rownames(mantel.comm.r)<-cor_colnames
  colnames(mantel.comm.r)<-cor_colnames
  mantel.comm.cil<-mantel.comm.ciu<-mantel.comm.p<-mantel.comm.r
  
  store_resultC<-function(i,j,MTR) {
    mantel.comm.cil[i,j] <<- MTR[5]
    mantel.comm.ciu[i,j] <<- MTR[6]
    mantel.comm.p[i,j] <<- MTR[4]
    mantel.comm.r[i,j] <<- MTR[1]
  }
  
  store_resultC(2,1,mantel(simdist ~ ibddist, nperm = 10000, mrank = T))
  store_resultC(3,1,mantel(simdist ~ cultdistsc, nperm = 10000, mrank = T))
  store_resultC(3,2,mantel(ibddist ~ cultdistsc, nperm = 10000, mrank = T))
  
  for(i in 1:rows) for(j in 1:cols) mantel.comm.r[i,j]<-mantel.comm.r[j,i]
  for(i in 1:rows) for(j in 1:cols) mantel.comm.p[i,j]<-mantel.comm.p[j,i]
  
# Save and remove
  resave(mantel.comm.cil, mantel.comm.ciu, mantel.comm.p, mantel.comm.r,
         file = "Results/manteltests.RData")
  remove(mantel.comm.cil, mantel.comm.ciu, mantel.comm.p, mantel.comm.r)
  
  
# Insect foraging
  cor_colnames<-c('Rare alleles', 
                  'Mean IBD length','Insect foraging')
  mantel.insect.r<-matrix(1,nrow=rows,ncol=cols)
  rownames(mantel.insect.r)<-cor_colnames
  colnames(mantel.insect.r)<-cor_colnames
  mantel.insect.cil<-mantel.insect.ciu<-mantel.insect.p<-mantel.insect.r
  
  store_resultC<-function(i,j,MTR) {
    mantel.insect.cil[i,j] <<- MTR[5]
    mantel.insect.ciu[i,j] <<- MTR[6]
    mantel.insect.p[i,j] <<- MTR[4]
    mantel.insect.r[i,j] <<- MTR[1]
  }
  
  store_resultC(2,1,mantel(simdist ~ ibddist, nperm = 10000, mrank = T))
  store_resultC(3,1,mantel(simdist ~ cultdistsi, nperm = 10000, mrank = T))
  store_resultC(3,2,mantel(ibddist ~ cultdistsi, nperm = 10000, mrank = T))
  
  for(i in 1:rows) for(j in 1:cols) mantel.insect.r[i,j]<-mantel.insect.r[j,i]
  for(i in 1:rows) for(j in 1:cols) mantel.insect.p[i,j]<-mantel.insect.p[j,i]

# Save and remove
  resave(mantel.insect.cil, mantel.insect.ciu, mantel.insect.p, mantel.insect.r,
         file = "Results/manteltests.RData")
  remove(mantel.insect.cil, mantel.insect.ciu, mantel.insect.p, mantel.insect.r)
  
# Foraging
  cor_colnames<-c('Rare alleles', 
                  'Mean IBD length','Foraging')
  mantel.foraging.r<-matrix(1,nrow=rows,ncol=cols)
  rownames(mantel.foraging.r)<-cor_colnames
  colnames(mantel.foraging.r)<-cor_colnames
  mantel.foraging.cil<-mantel.foraging.ciu<-mantel.foraging.p<-mantel.foraging.r
  
  store_resultC<-function(i,j,MTR) {
    mantel.foraging.cil[i,j] <<- MTR[5]
    mantel.foraging.cil[i,j] <<- MTR[6]
    mantel.foraging.p[i,j] <<- MTR[4]
    mantel.foraging.r[i,j] <<- MTR[1]
  }
  
  store_resultC(2,1,mantel(simdist ~ ibddist, nperm = 10000, mrank = T))
  store_resultC(3,1,mantel(simdist ~ cultdistsf, nperm = 10000, mrank = T))
  store_resultC(3,2,mantel(ibddist ~ cultdistsf, nperm = 10000, mrank = T))
  
  for(i in 1:rows) for(j in 1:cols) mantel.foraging.r[i,j]<-mantel.foraging.r[j,i]
  for(i in 1:rows) for(j in 1:cols) mantel.foraging.p[i,j]<-mantel.foraging.p[j,i]

# Save and remove
  resave(mantel.foraging.cil, mantel.foraging.ciu, mantel.foraging.p, mantel.foraging.r,
         file = "Results/manteltests.RData")
  remove(mantel.foraging.cil, mantel.foraging.ciu, mantel.foraging.p, mantel.foraging.r)

## CLEAR WORKSPACE..............................................................s
  rm(list= ls()[!(ls() %in% c('ddyads.raw','ddyads', "%!in%"))])

  detach("package:ecodist", unload=TRUE)

# ---- 
