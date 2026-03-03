#############################################################################################################
#############################################################################################################
## Functions for calculating and plotting isotopic niches (SEAc) and TNW
## Import these functions, then see accompanying script for analysis.
## 
## 3 functions: 
## (1) kapow2 -> adapts SIBER kapow function to plot species and community; 
## (2) hullarea -> imports internal SIBER function calculate area of group (i.e., species) ellipses 
## (3) "RINI" -> calculates Relative Individual Niche Index (RINI) at species level
## from Sheppherd et al. 2018 (Ecol Letters) and reports total niche width of population (TNW), 95% ellipse 
## area for each "group" (i.e. ispecies; SEAc), and species niche index (Niche.Index). RINI is
## the proportion of the population niche (TNW) occupied by each inividuals niche (Area.Ellipse). Scales from
## 0-1, where 0 is complete specialization and 1 is complete generalist.
##
## Written by: PJM, 03/16/2021
#############################################################################################################
#############################################################################################################

#rm(list=ls())

require(tidyverse)
require(dplyr)
require(SIBER)
require(spatstat)

#############################################################################################################
#############################################################################################################
## Function 1: Revamped kapow function for plotting
#############################################################################################################
#############################################################################################################

kapow2 = function(Y=data, do.plot = TRUE) {
  
  Y$group <- factor(Y$group)
  
  ellY <- siberKapow(Y, isoNames = c("iso1","iso2"), group = "group", pEll = .95)
  bdry <- as.data.frame(ellY$bdry[[1]])
  
  p <- ggplot(ellY$ell.coords, 
              mapping = aes_string(x = 'X1', y = 'X2', color = 'group', fill = 'group')) +
    geom_polygon(alpha = 0.1) + theme_classic() +  theme(legend.position = "none") +
    labs(x=expression(delta*{}^13*"C (\u2030)")) +
    labs(y=expression(delta*{}^15*"N (\u2030)")) +
    geom_polygon(data = bdry, mapping = aes_string(x =bdry$x, y =bdry$y), 
                 inherit.aes=F, size=2, colour="black", fill=NA)+
    geom_point(data=Y, mapping = aes_string(x ="iso1", y ="iso2", colour="group"))
  
  if (do.plot) print(p)
  
  return(p)
  
}

#############################################################################################################
#############################################################################################################
# Function 2: hullarea function to calculate ellipse areas
#############################################################################################################
#############################################################################################################

hullarea <- function (x,y) {
  ne <- length(x)
  harea <- abs (0.5 * ( (x[1:(ne-1)] %*% y[2:ne]) - ( y[1:(ne-1)] %*% x[2:ne]) ) )
  harea
}

#############################################################################################################
#############################################################################################################
# Function 3: RINI function to calculate niche metrics
#############################################################################################################
#############################################################################################################

rini <- function (z = data) { 
  
  # calculate total niche width -> area of all ellipses
  tnw <- hullarea(z$bdry[[1]]$x, z$bdry[[1]]$y)
  
  # calculate area of each "group" (i.e. species) ellipse
  seac = NULL
  
  for(i in 1:length(z$owin.coords)){
    
    #print(i)
    temp.x = z$owin.coords[[i]]$bdry[[1]]$x
    temp.y = z$owin.coords[[i]]$bdry[[1]]$y
    seac = rbind(seac, hullarea(temp.x, temp.y))
    
  }
  
  ## calculate specialization index for each individual 
  n.index = NULL
  
  for(i in 1:length(z$owin.coords)){
    
    #print(i)
    temp.x = z$owin.coords[[i]]$bdry[[1]]$x
    temp.y = z$owin.coords[[i]]$bdry[[1]]$y
    n.index = rbind(n.index, hullarea(temp.x, temp.y)/tnw)
    
  }
  
  #convert matrcies to dataframes
  id <- as.data.frame(unique(z$ell.coords$group))
  tnw <- as.data.frame(tnw)
  seac <- as.data.frame(seac)
  n.index <- as.data.frame(n.index)
  
  #combine as table
  table <- cbind(id, tnw, seac, n.index)
  
  #rename indices
  colnames(table) <- c("Species", "TNW", "SEAc", "Niche.Index")
  
  #output
  table
  
}

#############################################################################################################
#############################################################################################################
# el fin
#############################################################################################################
#############################################################################################################