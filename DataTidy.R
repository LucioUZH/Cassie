# ==============================================================================
# Cassandra Gunasekaram, AIM UZH, PhD - Chimp genetics project, 21.02.2023
# Data tidy and exploration
# Last update: 21.02.2023
# ==============================================================================

# ______________________________________________________________________________
# Workspace preparations and data input----

  rm(list=ls())
  setwd("C:/Users/cassi/OneDrive - Universität Zürich UZH/PhD/Genetics/Analysis/Repository")
  library(dplyr)
  library(ggpubr)
  
# Not in function
  '%!in%' <- function(x,y)!('%in%'(x,y))

# ----
# ______________________________________________________________________________
# Generate IBD data-----
  
  dibd <- select(read.csv("Data/d.ibd.csv"), -1)
  
# Calculate average of maximum IBD length between individuals per population dyad
  indibd <- dibd %>%
    group_by(Sample.y, Sample.x) %>%
    summarise(Site.x = Site.x, Site.y = Site.y,
              maxibd = max(Length..bp.))
  indibd$dyadID <- paste0(pmin(indibd$Site.x, indibd$Site.y), "_",
                          pmax(indibd$Site.x, indibd$Site.y))
  popibd <- indibd %>%
    group_by(dyadID) %>%
    summarise(avgmaxibd = mean(maxibd))
  
  remove(indibd)
  
# Create nr individuals per site and per population dyad
  ibdind <- select(dibd, 2,9)
  colnames(ibdind) <- c("Sample.y", "Site.y")
  ibdind <- rbind(ibdind, select(dibd, c(Sample.y, Site.y)))
  ibdind <- distinct(ibdind, Sample.y, .keep_all = T) #206 samples
  nrind <- as.data.frame(ibdind %>% group_by(Site.y) %>%
                           summarise(nrind = n()))
  
  dibd <- merge(dibd, nrind, by.x = "Site.x", by.y = "Site.y")
  dibd <- merge(dibd, nrind, by.x = "Site.y", by.y = "Site.y")
  dibd$totalind <- dibd$nrind.x + dibd$nrind.y
  
  remove(ibdind, nrind)
  
# Add dyadID
  dibd$dyadID <- paste0(pmin(dibd$Site.x, dibd$Site.y), "_",
                       pmax(dibd$Site.y, dibd$Site.x))

# Calculate mean ibd, max ibd and add maxibdcorrect, avg.maxibd
  
  dibd <- dibd %>% group_by(dyadID) %>%
    mutate(meanibd = mean(Length..bp.)) %>%
    mutate(maxibd = max(Length..bp.)) %>%
    ungroup() %>%
    distinct(dyadID, .keep_all = T)
  dibd$maxibdcorrect <- dibd$maxibd / dibd$totalind
  dibd <- dibd%>%
    select(c(21, 22:24, 20))
  dibd <- merge(dibd, popibd, by = "dyadID")
  
  write.csv(dibd, "Data/d.tidyibd.csv")
  remove(popibd, dibd)
 
# ----
# ______________________________________________________________________________
# Generate Private allele data-----
  
  load("Data/d.snps.RData")
  
# Calculate nr of snps and proportion sequenced in popdyads
  # Function
  nrdyadssnps <- function(data){

    x <- data.frame(pop1 = rep(colnames(data), each = ncol(data)))
    x$pop2 <- rep(colnames(data), times = ncol(data))
    x$dyadID <- paste0(pmin(x$pop1, x$pop2), "_",
                       pmax(x$pop1, x$pop2))
    x$ones <- NA
    x$total <- NA
    c = 0
    for (i in 1:ncol(data)) {
      for (j in 1:ncol(data)) {
        c = c+1
        x[c,4] <- length(which(data[,i] == 1 & data[,j] ==1))
        x[c,5] <- length(which(data[,i] %in% c(0,1) & data[,j] %in% c(0,1)))
        print(x[c,])
        
      }
            
    }
    x <- distinct(x, dyadID, .keep_all = T)
    return(x)
  }
    
    d.dyadsnps <- nrdyadssnps(d.snps)
    d.dyadsnps$sim <- d.dyadsnps$ones / d.dyadsnps$total

    hist(d.dyadsnps$sim)
    hist(filter(d.dyadsnps, pop1 != pop2)$sim)
    
# Creat bins according to position in list
    d.dyadsnps <- droplevels(filter(d.dyadsnps, pop1 != pop2))
    d.dyadsnps <- d.dyadssnps[order(d.dyadsnps$sim),]
    d.dyadsnps$index <- 1:nrow(d.dyadsnps)
    d.dyadsnps$simbin <- NA
    boundaries <- quantile(d.dyadsnps$index,na.rm = T,probs = seq(0,1, by = 0.1))
    i = 0
    while(i <10){
      i = i+1
      d.dyadsnps$simbin <- ifelse(d.dyadsnps$index >=boundaries[i] & d.dyadsnps$index <boundaries[i+1],
                            i, d.dyadsnps$simbin)
    }
    d.dyadsnps$simbin[which(d.dyadsnps$index == 561)] <-  10

# Export and remove
    write.csv(d.dyadsnps, "Data/d.tidysnps.csv")
    remove(d.snps, d.dyadsnps, i, boundaries)

# ----
# ______________________________________________________________________________
# Generate dyad data-----

    dsite <- select(read.csv("Data/d.site.csv"), -1)
    dsite <- droplevels(filter(dsite, pop %!in% c("Azagny", "Bia", "Loma")))
    
# dyads - popluations variables
    ddyads <- data.frame(t(combn(unique(dsite$pop), 2)))
    ddyads <- merge(ddyads, dsite, by.x = "X1", by.y = "pop")
    ddyads <- merge(ddyads, dsite, by.x = "X2", by.y = "pop")
    ddyads <- select(ddyads, c(2,1,3:5,7:11,13,14))
    
    colnames(ddyads) <- c("population1", "population2",
                          "lat1", "lon1", "region1", "habitat1", "obstime1",
                          "lat2", "lon2", "region2", 
                          "habitat2", "obstime2")
    remove(dsite)
 
# dyad variables
  ddyads$dyadID <- paste0(pmin(ddyads$population1, ddyads$population2),
                          "_",
                          pmax(ddyads$population1, ddyads$population2))
  ddyads$dyadregion <- ifelse(ddyads$region1 == ddyads$region2,
                              paste0(ddyads$region1, ddyads$region2),
                              "acrossregions")
  ddyads$sameregion <- ifelse(ddyads$region1 == ddyads$region2, 1, 0)
  ddyads$samehabitat <- ifelse(ddyads$habitat1 == ddyads$habitat2, 1, 0)
  ddyads$dyadhabitat <- paste0(pmin(ddyads$habitat1, ddyads$habitat2),
                               pmax(ddyads$habitat1, ddyads$habitat2))
  
  ddyads$distance <- NA
  for (i in (1:length(ddyads$population1))) {
    d <- as.numeric(geosphere::distm(c(ddyads$lon1[i], ddyads$lat1[i]), c(ddyads$lon2[i], ddyads$lat2[i])))
    ddyads$distance[i] <- d
  }
  
# Add IBD and rare allele data
  
  dibd <- select(read.csv("Data/d.tidyibd.csv"), -1)
  ddyads <- merge(ddyads, dibd, by = "dyadID", all.x = T)
  ddyads[is.na(ddyads) == T] <- 0
  ddyads$presenceibd <- ifelse(ddyads$meanibd != 0, 1, 0)
  remove(dibd)
  
  dsnps <- select(read.csv("Data/d.tidysnps.csv"), -1)
  ddyads <- merge(ddyads, dsnps, by = "dyadID")
  
  remove(dsnps)
    
# Add culturalknowledge
  #replicate lines to match number of shared traits per dyad, 
  #repeat = number of tools; times = number of dyads 
  
  dcult <- read.csv("Data/d.culture.csv", sep = ";")
  dbehaviour <- distinct(dcult, behaviour, .keep_all = T)
  
  ndyads <- as.numeric(nrow(ddyads))
  nbehaviour <- as.numeric(nrow(dbehaviour))
  
  ddyads <- ddyads[rep(1:nrow(ddyads),each=nbehaviour),]
  
  ddyads$behaviour <- rep(dbehaviour$behaviour[1:nbehaviour],
                          times = ndyads)
  ddyads$category <- rep(dbehaviour$category[1:nbehaviour], times = ndyads)
  ddyads$complexity <- rep(dbehaviour$complexity[1:nbehaviour], times = ndyads)
  
  remove(dbehaviour, nbehaviour, ndyads)
  
# Add presence in each population
  dpresence <- select(dcult, c(1,2,5))
  ddyads <- merge(ddyads, dpresence, by.x = c("population1", "behaviour"), 
                  by.y = c("site", "behaviour"))
  ddyads <- merge(ddyads, dpresence, by.x = c("population2", "behaviour"), 
                  by.y = c("site", "behaviour"))
  
  remove(dpresence, dcult)
  
# Order and rename columns
  colnames(ddyads)
  ddyads <- ddyads[c(3,1,4:32,2,33:36)]
  colnames(ddyads)[35:36] <- c("presence1", "presence2")
  
# Shared variable
  for ( i in (1:nrow(ddyads))){
    if (is.na(ddyads$presence1[i]) == T |  is.na(ddyads$presence2[i]) == T){
      ddyads$shared[i] <- NA
    } else if (is.na(ddyads$presence1[i]) == F &  is.na(ddyads$presence2[i]) == F & 
               ddyads$presence1[i] == 1 & ddyads$presence2[i] == 1){
      ddyads$shared[i] <- 1
    } else{
      ddyads$shared[i] <- 0
    }
  }   
  
  write.csv(ddyads,"Data/d.dyads.csv") 

# ----  