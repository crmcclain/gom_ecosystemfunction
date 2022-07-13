# title: "Woodfall Ecosystem Function"
# author: "Craig R. McClain"

####load packages####
  require(dplyr)
  require(ggplot2)
  require(vegan)
  require(tidyr)
  require(stringr)
  library(ggrepel)
  require(gridExtra)
  library(tibble)
  library(kableExtra)
  library(data.table)
  require (MuMIn)
  require(nlme)
  require(car)
  library(viridis)  
  library(candisc)
  library(emmeans)
  require(plyr) 
  require(Hmisc) #weighted variance



####figure theme for the paper####
  theme_craig <- function () { 
    theme_bw(base_size=12) %+replace% 
      theme(
        # change stuff here
        axis.line = element_line(colour = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        legend.position="none",
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0))}

  
  
####load and process data#####
  setwd("~/Dropbox/Return of the Woodfall/Wood Falls/WF Data")
  logbytaxa<- data.frame(read.csv("woodfall_logbytaxa.csv", header=TRUE, stringsAsFactors=FALSE))
  logsummary <-  data.frame(read.csv("woodfall_logs.csv", header=TRUE, stringsAsFactors=FALSE))
  NCOL <- ncol(logbytaxa)
  
  #removing stations and species with xero total abundance
  logbytaxa2 <- logbytaxa %>%
    replace(is.na(.), 0) %>% #replace all NAs (missing data cells) with xeros
    mutate(rsum = rowSums(.[2:NCOL])) %>% #remove species (rows) with sums of xeros
    filter(rsum>0) %>%
    select(-rsum)
  
  logbytaxa2 <- logbytaxa2[, colSums(logbytaxa2 != 0) > 0] #remove stations (cols) with sums of zero
  
  #time to transpose for vegan
  # first remember the names
  n <- logbytaxa2$ID
  
  # transpose all but the first column (name)
  logbytaxa3 <- as.data.frame(t(logbytaxa2[,-1]))
  colnames(logbytaxa3) <- n
  logbytaxa3$log_id <- factor(row.names(logbytaxa3))
  
  logbytaxa3 <- logbytaxa3[row.names( logbytaxa3) != "X.1", , drop = FALSE]
  
  ####overall woodfall diversity time###
  logdiversity <- logbytaxa3 %>%
    mutate(Abundance = rowSums(logbytaxa3[,1:141]),                 
           S = specnumber(logbytaxa3[,1:141]),
           H = diversity(logbytaxa3[,1:141],index="shannon"),
           Simp = diversity(logbytaxa3[,1:141],index="simpson"),
           log10Abundance=log10(Abundance)) %>%
    select(Abundance, log10Abundance, S, H, Simp, log_id)
  
  logdiversity$log_id<-str_replace_all(logdiversity$log_id, "[.]", "-")
  
  logdiversity <- left_join(logdiversity, logsummary, by=c("log_id"="Log.ID"))  
  
  logdiversity$log10Mass <- log10(logdiversity$wood.mass.final)
  
  ####ecosystem function metrics###
  logbytaxa_xylo <- logbytaxa3 %>% 
    select(contains("xylo")) %>%
    select(-"xylo-unid")
  
  logdiversity <- logdiversity %>%
    mutate(xylo_Abundance = rowSums(logbytaxa_xylo),                 
           xylo_S = specnumber(logbytaxa_xylo),
           xylo_H = diversity(logbytaxa_xylo,index="shannon"),
           xylo_Simp = diversity(logbytaxa_xylo,index="simpson"),
           xylo_log10Abundance=log10(xylo_Abundance),
           log10consumed = log10(wood.mass.final-weight_final),
           perc.consumed = ((wood.mass.final-weight_final)/wood.mass.final)*100,
           xylo3 = logbytaxa_xylo$`xylo-3`,
           xylo2 = logbytaxa_xylo$`xylo-2`,
           xylo1 = logbytaxa_xylo$`xylo-1`,
           xylo5 = logbytaxa_xylo$`xylo-5`,
           xylo11 = logbytaxa_xylo$`xylo-11`,
           xylo10 = logbytaxa_xylo$`xylo-10`,
           xylo12 = logbytaxa_xylo$`xylo-12`,
           xylo16 = logbytaxa_xylo$`Xylo-16`,
           xylo7 = logbytaxa_xylo$`xylo-7`,
           xylo15 = logbytaxa_xylo$`xylo-15`,
           xylo13 = logbytaxa_xylo$`xylo-13`,
           xylo9 = logbytaxa_xylo$`xylo-9`,
           xylo_Abundance_noxylo1 = xylo_Abundance - xylo1,
           log10xylo1 =log10(xylo1),
           log10_xyloabund_noxylo1 =log10(xylo_Abundance_noxylo1),
           wood.type = as.factor(wood.type)) 
  
  logdiversity$hardness <- revalue(logdiversity$wood.type, 
                                   c("Celtis laevigata"="3910", 
                                     "Pinus echinata"="3070", 
                                     "Pinus elliottii"="3380",
                                     "Quercus rubra"="5430",
                                     "Quercus virginiana"="12920",
                                     "Carya illinoiensis"="8100",
                                     "Magnolia grandiflora"="4540",     
                                     "Quercus alba"="5990",
                                     "Salix nigra"="1920",
                                     "Taxodium distichum"="2270",
                                     "Ulmus americana"="2830"))

  
  logdiversity$hardness<-as.numeric(as.character(logdiversity$hardness))
  
  logdiversity_drop <- logdiversity %>%
    filter(!log_id %in% c("L-G-31","X-1")) #dropped 31 because dropped by ROV on seafloor
  
  
####bivalve position data####
  setwd("~/Dropbox/Return of the Woodfall/Projects/Ecosystem Function/Code")
  logposition = data.frame(read.csv(file = "WF1_L-G-19_Quarter1_Coordinates.csv", header = T))
  
  logposition$Xylo_Group <- as.factor (logposition$Xylo_Group)
  
  
  xylo_pallette <- scale_color_manual(values = c("X. sp.1" = "#440154FF",
                                                 "A. brava" = "#39568CFF",
                                                 "S. sp.3" = "#238A8DFF" , 
                                                 "X. sp. 11" = "#FDE725FF"))
  xylo_pallette2 <- scale_fill_manual(values =  c("X. sp.1" = "#440154FF",
                                                  "A. brava" = "#39568CFF",
                                                  "S. sp.3" = "#238A8DFF" , 
                                                  "X. sp. 11" = "#FDE725FF"))
  
  logposition$Xylo_Group <- revalue(logposition$Xylo_Group, c("1"="X. sp.1", 
                                                              "2"="A. brava",
                                                              "3"= "S. sp.3",
                                                              "11" = "X. sp. 11"))
  
  ###manova
  xylo.man <- manova(cbind(Width_.mm., X_.mm., Y2, Z_.mm.) 
                     ~ Xylo_Group, data = logposition)
  summary.aov(xylo.man)

  
  ####emmeans plot#
  
  width_aov <- aov(Width_.mm.~Xylo_Group, data = logposition)
  width_aov.sum <- summary(width_aov)
  width_aov.emmeans <- emmeans(width_aov , ~ Xylo_Group)
  
  pw <- pwpp(width_aov.emmeans)+
    theme_craig()+
    xlab("P-Value")+
    ylab("Species")+
    xylo_pallette+
    xylo_pallette2+
    ggtitle("D. Burrow Width (mm)")+
    theme(plot.title = element_text(size = 10, face = "bold"))+
    annotate("rect", xmin = 0.05, xmax = 1.05, ymin = 0, ymax = 5.5, 
             alpha = .2)
  
  x_aov <- aov(X_.mm.~Xylo_Group, data = logposition)
  x_aov.sum <- summary(x_aov)
  x_aov.emmeans <- emmeans(x_aov , ~ Xylo_Group)
  
  px <- pwpp(x_aov.emmeans)+
    theme_craig()+
    xlab("P-Value")+
    ylab("Species")+
    xylo_pallette+
    xylo_pallette2+
    ggtitle("C. Depth From Surface, X distance (mm)")+
    theme(plot.title = element_text(size = 10, face = "bold"))+
    annotate("rect", xmin = 0.05, xmax = 1.05, ymin = 0, ymax = 5.5, 
             alpha = .2)
  
  
  y_aov <- aov(Y2~Xylo_Group, data = logposition)
  y_aov.sum <- summary(y_aov)
  y_aov.emmeans <- emmeans(y_aov , ~ Xylo_Group)
  
  py <- pwpp(y_aov.emmeans)+
    theme_craig()+
    xlab("P-Value")+
    ylab("Species")+
    xylo_pallette+
    xylo_pallette2+
    ggtitle("B. Proxmity to End, Y distance (mm)")+
    theme(plot.title = element_text(size = 10, face = "bold"))+
    annotate("rect", xmin = 0.05, xmax = 1.05, ymin = 0, ymax = 5.5, 
             alpha = .2)
  
  
  
  z_aov <- aov(Z_.mm.~Xylo_Group, data = logposition)
  z_aov.sum <- summary(z_aov)
  z_aov.emmeans <- emmeans(z_aov , ~ Xylo_Group)
  
  pz <- pwpp(z_aov.emmeans)+
    theme_craig()+
    xlab("P-Value")+
    ylab("Species")+
    xylo_pallette+
    xylo_pallette2+
    ggtitle("A. Surface Preference, Z distance (mm)")+
    theme(plot.title = element_text(size = 10, face = "bold"))+
    annotate("rect", xmin = 0.05, xmax = 1.05, ymin = 0, ymax = 5.5, 
             alpha = .2)
  
  pdf(file="figure6.pdf",width=8, height=6)
  par(mar=c(5,3,2,2)+0.1) #removes space from around edges of pdf
  grid.arrange(pz,py,px,pw)
  dev.off()   
  
  
####phylogenetic diversity####
  library(picante)
  phy <- read.nexus("XyloTeredi.nex")
  class(phy)
  phy$tip.label <- c("xylo-7", "Teredo_clappi", 
                     "xylo-15" ,"xylo-13",
                     "xylo-10", "xylo-12", 
                     "xylo-16", "xylo-11" ,
                     "xylo-5", "xylo-2", 
                     "xylo-3", "xylo-1",
                     "xylo-4", "xylo-1A", 
                     "xylo-1B" , "xylo-1C" ,
                     "xylo-1D")
  Ntip(phy)
  plot(phy, cex = 0.5)
  
  # check for mismatches/missing species
  combined <- match.phylo.comm(phy, logbytaxa_xylo)
  # the resulting object is a list with $phy and $comm elements.  replace our
  # original data with the sorted/matched data
  phy <- combined$phy
  comm <- combined$comm
  
  # Calculate Faith's PD
  comm.pd <- pd(comm, phy)
  comm.pd$wood_id <- row.names(comm.pd)
  comm.pd$wood_id <-gsub("[.]", "-", comm.pd$wood_id)
  
  
  
  logdiversity_drop <- left_join(logdiversity_drop, comm.pd, by=c("log_id"="wood_id")) %>%
    select(Abundance, S, H, Simp, log_id, site,log10Mass:hardness, PD, wood.type)
  
  glimpse(logdiversity_drop)
  
  
####Species Summary Table####
  
  logposition.sum <- logposition %>% 
    dplyr::group_by(Xylo_Group) %>%
    dplyr::summarise(Surface.Preference=mean(Z_.mm.),
                     Proxmity.to.End=mean(Y_.mm.),
                     Depth.from.Surface=mean(X_.mm.),
                     Burrow.Width=mean(Width_.mm.)) %>%
    arrange(desc(Burrow.Width))
  
  logposition.sum$Xylo_Group <- as.factor(logposition.sum$Xylo_Group)
  
  
  xylo_summary <- logbytaxa_xylo %>%
    replace(is.na(.), 0) %>%
    summarise_all(funs(sum)) 
  
  xylo_summary[2, 1:12] <- colSums(logbytaxa_xylo != 0)
  xylo_summary <- as_tibble(t(xylo_summary), rownames = "species")
  colnames(xylo_summary) <- c("Species", "Abundance", "Number of Wood Falls")
  
  xylo_summary$Species_name<- revalue(xylo_summary$Species, c("xylo-3"="Spiniapex sp. 3",
                                                              "xylo-2"="Abidtoconus cf. brava",
                                                              "xylo-1"="Xyloredo sp.1",
                                                              "xylo-5"="Xylopholas sp. 5",
                                                              "xylo-11"="Xylophagidae sp. 11",
                                                              "xylo-10"="Xylonora sp.10",
                                                              "xylo-12"="Xylonora sp. 12",
                                                              "Xylo-16"="Xylophagidae sp. 16",
                                                              "xylo-7"="Bankia sp. 7",
                                                              "xylo-15"="Xylonora sp. AC78",
                                                              "xylo-13"="Xylonora sp. 13",
                                                              "xylo-9" ="Xylophaga sp. 9"))
  
  
  
  xylo_summary <- xylo_summary%>%
    arrange(desc(Abundance))
  
  xylo_summary$Size_mg <- c(125.93,28.01 ,1.13,28.48,67.42,12.88,18.44,0.55,2.16,151.67,82.82,15.90)   
  
  logposition.sum %>%
    kable(caption = "Table: Summary Statistics for Xylophagid Bivalve Species. All measurements in mm", digits=2, col.names = c("Xylophagid Species", "Surface Preference", "Proximity to End", "Depth from Surface", "Burrow Width"))%>%
    kable_classic_2(full_width = F)
  

####Create Body Size Analysis####
  
  #create a triplet
  
  logbytaxa_xylo$wood_id <- row.names(logbytaxa_xylo)
  
  xylo_triplet <- logbytaxa_xylo %>%
    gather(species, abundance, -wood_id) %>%
    filter(abundance > 0)
  
  xylo_triplet <- left_join(xylo_triplet, xylo_summary, by = c("species" = "Species")) 
  
  #get summary metabolic rate per station
  xylo_triplet <- xylo_triplet %>%
    mutate(biomass=abundance*Size_mg) %>%
    select(-Abundance, -`Number of Wood Falls`)
  
  
  
  size_summary <- xylo_triplet  %>%
    dplyr::group_by(wood_id) %>%
    dplyr::summarise(tot_biomass=sum(biomass),
                     mean_sz = mean(log10(Size_mg)),
                     max_sz = max(log10(Size_mg)),
                     min_sz = min (log10(Size_mg)),
                     var_sz = var(log10(Size_mg)),
                     w_mean_sz = weighted.mean(log10(Size_mg), abundance),
                     w_var_sz = wtd.var(log10(Size_mg), abundance),
                     sd_sz = sd(log10(Size_mg)),
                     w_sd_size =sqrt(w_var_sz))
  
  size_summary$wood_id <-gsub("[.]", "-", size_summary$wood_id )
  
  logdiversity_drop <- left_join(logdiversity_drop, size_summary, by=c("log_id"="wood_id"))  
  

  
  
####Figure 1 Species Summary Plot####
  p1<- ggplot(logposition, aes(group=Xylo_Group, x=Z_.mm., fill=Xylo_Group, color=Xylo_Group)) +
    geom_density(alpha=0.5) + 
    theme_craig()+
    ylab("Density")+
    xlab("Surface Preference, Z distance (mm)")+
    xylo_pallette+
    xylo_pallette2+
    ggtitle("A")
  
  p2<-ggplot(logposition, aes(group=Xylo_Group, x=Y2, fill=Xylo_Group, color=Xylo_Group)) +
    geom_density(alpha=0.5) + 
    theme_craig()+
    ylab("Density")+
    xlab("Proxmity to End, Y distance (mm)")+
    xylo_pallette+
    xylo_pallette2+
    ggtitle("B")
  
  p3<-ggplot(logposition, aes(group=Xylo_Group, x=X_.mm., fill=Xylo_Group, color=Xylo_Group)) +
    geom_density(alpha=0.5) + 
    theme_craig()+
    ylab("Density")+
    xlab("Depth From Surface, X distance (mm)")+
    xylo_pallette+
    xylo_pallette2+
    ggtitle("C")
  
  p4<-ggplot(logposition, aes(group=Xylo_Group, x=Width_.mm., fill=Xylo_Group, color=Xylo_Group)) +
    geom_density(alpha=0.5) + 
    theme_craig()+
    ylab("Density")+
    xlab("Burrow Width (mm)")+
    xylo_pallette+
    xylo_pallette2+
    ggtitle("D")
  
  
  xylo_pallette3 <- scale_color_manual(values = c("Xyloredo sp.1" = "#440154FF",
                                                 "Abidtoconus cf. brava" = "#39568CFF",
                                                 "Spiniapex sp. 3" = "#238A8DFF" , 
                                                 "Xylophagidae sp. 11" = "#FDE725FF"))
  
  xylo_pallette4 <- scale_fill_manual(values =  c("Xyloredo sp.1" = "#440154FF",
                                                  "Abidtoconus cf. brava" = "#39568CFF",
                                                  "Spiniapex sp. 3" = "#238A8DFF" , 
 
                                                                                                   "Xylophagidae sp. 11" = "#FDE725FF"))
 
  
  p5<-ggplot(xylo_triplet, aes(x=(reorder(Species_name,abundance)), y=abundance, fill=Species_name, color=Species_name)) +
    geom_boxplot(alpha=0.5) + 
    geom_jitter(width=0.1, alpha=0.5, pch=21, cex=3, color="grey20")+
    theme_craig()+
    scale_y_log10()+
    xlab("Species")+
    ylab("Abundance")+
    xylo_pallette3+
    xylo_pallette4+
    theme(axis.text.x = element_text(angle = 45, h=1))
    
 
  library(ggpubr)

  
  pdf(file="figure1.pdf",width=8, height=6)
  par(mar=c(5,3,2,2)+0.1) #removes space from around edges of pdf
  ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
  dev.off()   
  
  
  pdf(file="figure3.pdf",width=8, height=6)
  par(mar=c(5,3,2,2)+0.1) #removes space from around edges of pdf
  p5
  dev.off()   
  
  
####Abundance AOV####
  abundance_model <- aov(log10(abundance)~Species_name, data=xylo_triplet)
  summary(abundance_model)
  abundance_model
  
  abundance_aov.emmeans <- emmeans(abundance_model , ~ Species_name)
  
  a1 <- pwpp(abundance_aov.emmeans)+
    theme_craig()+
    xlab("P-Value")+
    ylab("Species")+
    xylo_pallette3+
    xylo_pallette4+
    ggtitle("Log10 Abundnace")+
    theme(plot.title = element_text(size = 10, face = "bold"))+
    annotate("rect", xmin = 0.05, xmax = 1.05, ymin = 0, ymax =13, 
             alpha = .2)
  
  a1
  
####species null models####
  
  #creating dropped species by log 
  logbytaxa_xylo_drop <- logbytaxa_xylo %>%
    mutate(wood_id=row.names.data.frame(logbytaxa_xylo)) %>%
    filter(!wood_id %in% c("L.G.31","X.1", "L.Y.29")) %>%
    select(-wood_id)
  #permutation size and number of species
  Species_no <- ncol(logbytaxa_xylo_drop)
  permutations <- 1000
  
  
  #create species summary data frame to hold
  Species_sum<-data.frame(matrix(NA, nrow = Species_no, ncol = 1))
  colnames(Species_sum)[1]<-"Species.Name"
  Species_sum$Species.Name <- colnames(logbytaxa_xylo_drop)
  glimpse(Species_sum)
  
  #create data frame of permuted function (wood consumed) data
  VarPerm <- data.frame(replicate(permutations, slice(data.frame(logdiversity_drop$log10consumed), sample(1:n()))))
  
  for (xylo_loop in 1:Species_no) {
    
    
    #calculate variables for presence/absence variables for permuted data
    
    VarPerm2<-data.frame(matrix(NA, nrow = permutations, ncol = 1))
    
    VarPerm2$SoA <- colMeans(VarPerm[(logbytaxa_xylo_drop[,xylo_loop]==0),1:permutations],na.rm = TRUE)
    VarPerm2$SoP <- colMeans(VarPerm[(logbytaxa_xylo_drop[,xylo_loop]>0),1:permutations],na.rm = TRUE)
    VarPerm2$Diff <- VarPerm2$SoP-VarPerm2$SoA
    
    
    #calculate regression data for abundance permuted data
    for (i in 1:permutations) {
      reg<- lm(logbytaxa_xylo_drop[,xylo_loop] ~ VarPerm[,i])
      VarPerm2$slope[i] <- reg$coefficients[2]
    }
    
    Species_temp <- NULL
    
    Species_temp <- matrix(NA, nrow = 1, ncol = 1)
    #finish off final metrics for species table
    
    #empirical for species presence/absence
    Species_temp$SoA       = mean(logdiversity_drop$log10consumed[(logbytaxa_xylo_drop[,xylo_loop]==0)], na.rm = TRUE)
    Species_temp$SoP       = mean(logdiversity_drop$log10consumed[(logbytaxa_xylo_drop[,xylo_loop]>0)], na.rm = TRUE)
    Species_temp$Diff      = Species_temp$SoP-Species_temp$SoA
    
    #empirical for species regression
    Species_temp$Slope     = lm(logdiversity_drop$log10consumed~logbytaxa_xylo_drop[,xylo_loop])$coefficients[2]
    
    
    #simulated for species presence/absence
    Species_temp$Diff_mean_sim   = mean(VarPerm2$Diff,na.rm = TRUE)
    Species_temp$Diff_mean_stdev = sd(VarPerm2$Diff,na.rm = TRUE)
    Species_temp$Diff_z_value    = (Species_temp$Diff-Species_temp$Diff_mean_sim)/Species_temp$Diff_mean_stdev
    Species_temp$Diff_p_value    = sprintf("%.04f",(min(c(sum(Species_temp$Diff  > VarPerm2$Diff),sum(Species_temp$Diff  < VarPerm2$Diff)))/permutations))
    
    #simulated for species regression
    Species_temp$Slope_mean_sim  = mean(VarPerm2$slope,na.rm = TRUE)
    Species_temp$Slope_mean_stdev  = sd(VarPerm2$slope,na.rm = TRUE)
    Species_temp$Slope_z_value    = (Species_temp$Slope-Species_temp$Slope_mean_sim)/Species_temp$Slope_mean_stdev
    Species_temp$Slope_p_value    = sprintf("%.04f",(min(c(sum(Species_temp$Slope  > VarPerm2$slope),sum(Species_temp$Slope  < VarPerm2$slope)))/permutations))
    
    
    Species_sum[xylo_loop,2:13] <- data.frame(Species_temp[2:13]) 
    
  }
  
  
  Species_sum <- Species_sum %>%
    select(-SoA, -SoP, -Diff_mean_stdev, -Slope_mean_stdev) %>%
    relocate (Slope, .after = Diff_p_value) 
  
  
  xylo_summary[,6:13] <- Species_sum[,2:9]
  
  
  
  xylo_summary%>% 
    select(-Species) %>%
    relocate (Species_name, .before=Abundance) %>%
    kbl(caption="Table: Summary Statistics for Xylophagid Bivalve Species. 
           Randomization test for quantifying species importance to ecosystem function", 
        col.name = (c("Species Name", 
                      "Abundance", "Number of Wood Falls", "Mass (gm)",
                      "Empirical", "Simulated" , "Z-Score", "p-value",
                      "Empirical", "Simulated" , "Z-Score", "p-value")),
        digits = c(0, 0, 0, 1, 3, 3, 2, 4, 3, 3, 2, 4)) %>%
    kable_classic(full_width = F, html_font = "Cambria") %>%
    add_header_above(c(" " = 4, "Presence/Absence" = 4, "Slope" = 4)) %>%
    column_spec(7, color = ifelse(xylo_summary$Diff_z_value > 2, "blue", 
                                  ifelse(xylo_summary$Diff_z_value < -2, "red", "grey")))%>%
    column_spec(8, color = ifelse(xylo_summary$Diff_z_value > 2, "blue", 
                                  ifelse(xylo_summary$Diff_z_value < -2, "red", "grey")))%>%
    column_spec(11, color = ifelse(xylo_summary$Slope_z_value  > 2, "blue", 
                                   ifelse(xylo_summary$Slope_z_value < -2, "red", "grey")))%>%
    column_spec(12, color = ifelse(xylo_summary$Slope_z_value  > 2, "blue", 
                                   ifelse(xylo_summary$Slope_z_value  < -2, "red", "grey")))
  
  
  cor.test(log10(xylo_summary$Abundance), xylo_summary$Diff_z_value)
  cor.test(log10(xylo_summary$`Number of Wood Falls`), xylo_summary$Diff_z_value)
  

  
  
#####Figure 3 #####
  xylo_A <- ggplot(data=logdiversity_drop, aes(x=xylo_Abundance, log10consumed))+
    geom_point(cex=4, pch=21, alpha=0.5, col="black", fill="grey70")+
    geom_smooth(method=lm, se=FALSE, color="grey70")+
    #geom_text(aes(label=log_id))+
    theme_craig()+
    ylab("Consumed Weight (kg)")+
    xlab("Abundance")+
    ggtitle("A")
  
  xylo_S <- ggplot(data=logdiversity_drop, aes(x=xylo_S, y=log10consumed))+
    geom_point(cex=4, pch=21, alpha=0.5, col="black", fill="grey70")+
    geom_smooth(method=lm, se=FALSE, color="grey70")+
    theme_craig()+
    #geom_text(aes(label=log_id))+
    ylab("Consumed Weight (kg)")+
    scale_x_continuous(breaks=c(0,2,4,6,8,10, 12))+
    xlab("S")+
    ggtitle("B")
  
  xylo_J <- ggplot(data=  logdiversity_drop, aes(x=xylo_Simp,y=log10consumed))+
    geom_point(cex=4, pch=21, alpha=0.5, col="black", fill="grey70")+
    theme_craig()+
    ylab("Consumed Weight (kg)")+
    xlab("Simpson's J")+
    ggtitle("C")
  
  xylo_PD <- ggplot(data=logdiversity_drop, aes(x=PD,log10consumed))+
    geom_point(cex=4, pch=21, alpha=0.5, col="black", fill="grey70")+
    theme_craig()+
    geom_smooth(method=lm, se=FALSE, color="grey70")+
    ylab("Consumed Weight (kg)")+
    xlab("Faith's PD")+
    ggtitle("D")
  
  xylo_sd <- ggplot(data=logdiversity_drop, aes(x=sd_sz, log10consumed))+
    geom_point(cex=4, pch=21, alpha=0.5, col="black", fill="grey70")+
    geom_smooth(method=lm, se=FALSE, color="grey70")+
    #geom_text(aes(label=log_id))+
    theme_craig()+
    ylab("Consumed Weight (kg)")+
    xlab("Standard Deviation in Size")+
    ggtitle("E")
  
  xylo_mean_sz <- ggplot(data=logdiversity_drop, aes(x=mean_sz,log10consumed))+
    geom_point(cex=4, pch=21, alpha=0.5, col="black", fill="grey70")+
    theme_craig()+
    #geom_text(aes(label=log_id))+
    geom_smooth(method=lm, se=FALSE, color="grey70")+
    ylab("Consumed Weight (kg)")+
    xlab("Mean Size")+
    ggtitle("F")
  
  
  setwd("~/Dropbox/Return of the Woodfall/Projects/Ecosystem Function/New Paper")
  pdf(file="figure3_updated.pdf",width=11.5, height=8)
  par(mar=c(5,3,2,2)+0.1) #removes space from around edges of pdf
  grid.arrange(xylo_A, xylo_S, xylo_J, xylo_PD,xylo_sd, xylo_mean_sz)
  dev.off()    

  

  


  
###ecosystem function models####
  require(performance)
  #abundance

  glm_A_hard_interaction <- glm(log10consumed ~ xylo_Abundance+hardness+xylo_Abundance*hardness ,data=logdiversity_drop)
    summary(glm_A_hard_interaction)
    glm_A_hard<- glm(log10consumed ~ xylo_Abundance+hardness,data=logdiversity_drop)
    summary(glm_A_hard)
    
    #abundance reduced model and comparisons
    glm_A <- glm(log10consumed ~ xylo_Abundance,data=logdiversity_drop)
    summary(glm_A)
    AIC(glm_A_hard_interaction, glm_A_hard,glm_A)
    
  #richness
  
  glm_S_hard_interaction <- glm(log10consumed ~ xylo_S+hardness+ xylo_S*hardness,data=logdiversity_drop)
    summary(glm_S_hard_interaction)
    glm_S_hard <- glm(log10consumed ~ xylo_S+hardness,data=logdiversity_drop)
    summary(glm_S_hard)
  
    #richness reduced model and comparisons  
    glm_S <- glm(log10consumed ~ xylo_S,data=logdiversity_drop)
    summary(glm_S)
    AIC(glm_S_hard_interaction,  glm_S_hard, glm_S)
    
   #PD
    
    glm_PD_hard_interaction <- glm(log10consumed ~ PD+hardness+ PD*hardness,data=logdiversity_drop)
    summary(glm_PD_hard_interaction)
    glm_PD_hard <- glm(log10consumed ~ PD+hardness,data=logdiversity_drop)
    summary(glm_PD_hard)
    
    glm_PD <- glm(log10consumed ~ PD,data=logdiversity_drop)
    summary(glm_PD)


  #simpson's J and H
  summary(glm(log10consumed ~ xylo_H,data=logdiversity_drop))
  summary(glm(log10consumed ~ xylo_Simp,data=logdiversity_drop))
  
  
  
#####overall comparison of models####


 glm_size <- glm(log10consumed ~ sd_sz,data=logdiversity_drop)
 glm_S_A_PD_size <- glm(log10consumed ~ xylo_Abundance + xylo_S + PD + sd_sz,data=logdiversity_drop)
 glm_S_A_hard_PD_sz <-glm(log10consumed ~ xylo_Abundance + xylo_S + xylo_S*hardness + sd_sz +PD,data=logdiversity_drop)
 summary(glm_S_A_PD_size )
 
 glm_S_A_hard_sz <-glm(log10consumed ~ xylo_Abundance + xylo_S + xylo_S*hardness + sd_sz,data=logdiversity_drop)

 glm_A_size <- glm(log10consumed ~ sd_sz + xylo_Abundance,data=logdiversity_drop)
 

 #model table
 
     model_table <- AIC(glm_A, 
                        glm_S, 
                        glm_PD, 
                        glm_size,
                        glm_S_A_PD_size,
                        glm_S_A_hard_PD_sz,
                        glm_S_A_hard_sz,
                        glm_A_size)
    
     model_table$R2[1] <- r2_nagelkerke(glm_A)
     model_table$R2[2] <- r2_nagelkerke(glm_S)
     model_table$R2[3] <- r2_nagelkerke(glm_PD)
     model_table$R2[4] <- r2_nagelkerke(glm_size)
     model_table$R2[5] <- r2_nagelkerke(glm_S_A_PD_size)
     model_table$R2[6] <- r2_nagelkerke(glm_S_A_hard_PD_sz)
     model_table$R2[7] <- r2_nagelkerke(glm_S_A_hard_sz)
     model_table$R2[8] <- r2_nagelkerke(glm_A_size)
    
    
    row.names(model_table) <- c("Abundance Only", "Richness Only", "Phylogenetic Only", "Size Only", 
                                "Abundance, Richness, Phylogenetic, Size",
                                "Abundance, Richness, Phylogenetic, Size, Wood Hardness", 
                                "Abundance, Richness, Size, Wood Hardness",
                                "Abundance and Size")
     
     model_table  %>%
       kable(caption = "Table: Summary statistics for models", digits=2) %>%
       kable_classic_2(full_width = F)
     
     #effect sizes
     
     require(effectsize)
  
     model_sum <- data.frame(summary(glm_A_size)$coefficients)
     ETA <- data.frame(epsilon_squared(aov(glm_A_size))[2])
     model_sum$Epsilon2_partial[1] <- NA
     model_sum$Epsilon2_partial[2:3] <- ETA$Epsilon2_partial
     
     #row.names(model_sum) <- c("Intercept", "Abundance", "Richness","Phylogenetic Diveristy", "St. Dev. Body Size")
     
     model_sum%>%
       kable(caption = "Table: Summary statistics for model", digits =4) %>%
       kable_classic_2(full_width = F)
     

    
#####basic stats for paper#####
    summary(lm(log10consumed~xylo_Abundance, data = logdiversity_drop))
    summary(lm(log10consumed~xylo_S, data = logdiversity_drop))
    summary(lm(log10consumed~xylo_Simp, data = logdiversity_drop))
    summary(lm(log10consumed~PD, data = logdiversity_drop))
    
    summary(lm(log10consumed~sd_sz, data = logdiversity_drop))
    summary(lm(log10consumed~mean_sz, data = logdiversity_drop))
  
    
    cor.stats <- cor.test(logdiversity_drop$xylo_S, logdiversity_drop$xylo_Abundance, method=c("pearson"))
    cor.stats
    
    glm_S_A_hard_PD_sz <-glm(log10consumed ~ xylo_Abundance + xylo_S + xylo_S*hardness + mean_sz + PD + site,data=logdiversity_drop)
    
    

####10 logs analysis####
    tenlogs <- logdiversity_drop %>%
      filter(log_id %in% c("L-Y-29",
                      "L-W-22",
                      "L-W-25",
                      "L-Y-31",
                      "L-W-28",
                      "L-W-34",
                      "L-G-16",
                      "L-G-10",
                      "L-G-7",
                      "L-G-34"))
    tenlogs$TNW <- c(12.75,
                     6.43,
                     26.43,
                     19.23,
                     18.7,
                     20.31,
                     18.91,
                     12.55,
                     17.75,
                     12.49)
    
    tenlogs$RNI <- c(1,
                     0.705,
                     0.805,
                     1,
                     0.606666667,
                     0.503333333,
                     0.4125,
                     0.4075,
                     0.418,
                     0.503333333)
    
    
    
    glm_S_10 <- glm(log10consumed ~ xylo_S,data=tenlogs)
    glm_A_10 <- glm(log10consumed ~ xylo_Abundance,data=tenlogs)
    glm_PD_10 <- glm(log10consumed ~ PD,data=tenlogs)
    glm_sz_10 <- glm(log10consumed ~ sd_sz,data=tenlogs)
    glm_TNW_10 <- glm(log10consumed ~ TNW,data=tenlogs)
    glm_RNI_10 <- glm(log10consumed ~ RNI,data=tenlogs)

    model_table <- AIC(glm_S_10, 
                       glm_A_10, 
                       glm_PD_10, 
                       glm_sz_10,
                       glm_TNW_10,
                       glm_RNI_10)
    
    model_table$R2[1] <- r2_nagelkerke(glm_S_10)
    model_table$R2[2] <- r2_nagelkerke(glm_A_10)
    model_table$R2[3] <- r2_nagelkerke(glm_PD_10)
    model_table$R2[4] <- r2_nagelkerke(glm_sz_10)
    model_table$R2[5] <- r2_nagelkerke(glm_TNW_10)
    model_table$R2[6] <- r2_nagelkerke(glm_RNI_10)

    
    
    row.names(model_table) <- c("Richness", "Abundance", "Phylogenetic", "St. Dev. Body Size", "Total Isotopic Niche Width", "Relative Niche Index" )
    
    model_table  %>%
      kable(caption = "Table: Summary statistics for models", digits=2) %>%
      kable_classic_2(full_width = F)
    

    cor.test(tenlogs$sd_sz, tenlogs$TNW, method=c("pearson"))
    cor.test(tenlogs$xylo_S, tenlogs$TNW, method=c("pearson"))
    cor.test(tenlogs$sd_sz, tenlogs$xylo_S, method=c("pearson"))
    
    
    

