# load packages and data
# estimate individual functions of spp from presence absences
library(tidyverse)
devtools::install_github("ctkremer/priceTools")
library(priceTools)
library(grid)
library(gridExtra)
library(brms)
library(mgcv)
rstan::rstan_options(auto_write = TRUE)
theme_set(theme_minimal())
# should all the pairwise price paritions be rerun? This may take some time
# if FALSE previously saved objects are loaded.
rerun = FALSE

####!!!-- HELPER FUNCTIONS --!!!####

set_site_order = function(x,...){
  max_S = which(x$xylo_S == max(x$xylo_S))
  max_f = which(x$log10consumed == max(x$log10consumed))
  max_n = which(x$xylo_Abundance == max(x$xylo_Abundance))
  
  if(!identical(unlist(x$log10consumed)[1],unlist(x$log10consumed)[2])){
    x$siteNumber = 'site2'
    x[max_f,'siteNumber'] = 'site1'
  } else if(!identical(unlist(x$max_n)[1],unlist(x$max_n)[2])){
    x$siteNumber = 'site2'
    x[max_n, 'siteNumber'] = 'site1'
  } else{
    x$siteNumber = 'site2'
    x[max_S, 'siteNumber'] = 'site1'
  }
  x
}

set_site_order_n = function(x,...){
  max_S = which(x$xylo_S == max(x$xylo_S))
  max_f = which(x$log10consumed == max(x$log10consumed))
  max_n = which(x$xylo_Abundance == max(x$xylo_Abundance))
  
  if(!identical(unlist(x$xylo_Abundance)[1],unlist(x$xylo_Abundance)[2])){
    x$siteNumber = 'site2'
    x[max_n,'siteNumber'] = 'site1'
  } else if(!identical(unlist(x$log10consumed)[1],unlist(x$log10consumed)[2])){
    x$siteNumber = 'site2'
    x[max_f, 'siteNumber'] = 'site1'
  } else{
    x$siteNumber = 'site2'
    x[max_S, 'siteNumber'] = 'site1'
  }
  x
}

make_price_list = function(x,...){
  siteIndex = which(names(logList) %in% x)
  priceList = list(logList[[siteIndex[1]]], logList[[siteIndex[2]]])
  priceList = priceList %>% 
    purrr::map(~.x %>% dplyr::select(species, func = "consumption"))
  priceList
}

summarise_price = function(pricePart,...){
  priceComps = pricePart[[1]]
  deltaF = abs(priceComps["x.func"] - priceComps["y.func"])
  deltaS = abs(priceComps["x.rich"] - priceComps["y.rich"])
  
  data.frame(deltaF = deltaF,
             deltaS = deltaS,
             sharedS = priceComps["c.rich"],
             SRE.L = priceComps["SRE.L"],
             SRE.G = priceComps["SRE.G"],
             SIE.L = priceComps["SIE.L"],
             SIE.G = priceComps["SIE.G"],
             CDE = priceComps["CDE"],
             SL = priceComps["SL"],
             SG = priceComps["SG"],
             SR = priceComps["SR"],
             CE = priceComps["CE"]) %>%
    dplyr::mutate(across(matches("^\\w{2}$"), list(rel = ~(.x-deltaF)/deltaF)))
  
}
# internal from priceTools
process.data.price_cor = function(data, group.vars = NULL, standardize = TRUE){
  if(standardize){
    comps <- c("SRE.L", "SRE.G", "SIE.L", 
               "SIE.G", "CDE")
    data[, comps] <- 100 * (data[, comps]/data$x.func)
    dataStandard <- data
    data$y.func <- 100 * (data$y.func - data$x.func)/data$x.func
    data$x.func <- 0
    data$SIE.L <- data$SRE.L + data$SIE.L
    data$SRE.G <- data$SIE.L + data$SRE.G
    data$SIE.G <- data$SRE.G + data$SIE.G
  }
  else {
    data$x.func <- data$x.func
    data$y.func <- data$y.func
    data$SRE.L <- data$x.func + data$SRE.L
    data$SIE.L <- data$SRE.L + data$SIE.L
    data$SRE.G <- data$SIE.L + data$SRE.G
    data$SIE.G <- data$SRE.G + data$SIE.G
  }
  cols <- c(group.vars, "x.func", "SRE.L", "SIE.L", 
            "SRE.G", "SIE.G", "y.func", "x.rich", 
            "c.rich", "y.rich")
  p2 <- reshape2::melt(data[, cols], id.vars = c(group.vars, 
                                                 "x.rich", "c.rich", "y.rich"))
  p2$rich <- ifelse(p2$variable == "x.func", p2$x.rich, 
                    ifelse(p2$variable == "SRE.L", p2$c.rich, ifelse(p2$variable == 
                                                                       "SIE.L", p2$c.rich, ifelse(p2$variable == "SRE.G", 
                                                                                                  p2$y.rich, ifelse(p2$variable == "SIE.G", p2$y.rich, 
                                                                                                                    ifelse(p2$variable == "y.func", p2$y.rich, 
                                                                                                                           NA))))))
  if (!is.null(group.vars)) {
    p3b <- p2 %>% group_by_(.dots = c(group.vars, "variable")) %>% 
      summarise(mean.y = mean(value), y.qt.lw = quantile(value, 
                                                         probs = 0.025), y.qt.up = quantile(value, probs = 0.975), 
                mean.x = mean(rich), x.qt.lw = quantile(rich, 
                                                        probs = 0.025), x.qt.up = quantile(rich, probs = 0.975))
  }  else {
    p3b <- p2 %>% group_by(variable) %>% summarise(mean.y = mean(value), 
                                                   y.qt.lw = quantile(value, probs = 0.025), y.qt.up = quantile(value, 
                                                                                                                probs = 0.975), mean.x = mean(rich), x.qt.lw = quantile(rich, 
                                                                                                                                                                        probs = 0.025), x.qt.up = quantile(rich, probs = 0.975))
  }
  p4 <- p3b
  p2$variable <- factor(p2$variable, levels = c("x.func", 
                                                "SRE.L", "SIE.L", "SRE.G", "SIE.G", 
                                                "y.func"), labels = c("baseline", "SRE.L", 
                                                                      "SIE.L", "SRE.G", "SIE.G", "comparison"))
  p3b$variable <- factor(p3b$variable, levels = c("x.func", 
                                                  "SRE.L", "SIE.L", "SRE.G", "SIE.G", 
                                                  "y.func"), labels = c("baseline", "SRE.L", 
                                                                        "SIE.L", "SRE.G", "SIE.G", "comparison"))
  p4$variable <- as.character(p4$variable)
  p4$variable <- ifelse(p4$variable == "y.func", "SIE.G", 
                        p4$variable)
  p4$variable <- factor(p4$variable, levels = c("x.func", 
                                                "SRE.L", "SIE.L", "SRE.G", "SIE.G"), 
                        labels = c("SRE.L vector", "SIE.L vector", 
                                   "SRE.G vector", "SIE.G vector", "CDE vector"))
  p2$variable <- factor(p2$variable, levels = c("baseline", 
                                                "SRE.L", "SIE.L", "SRE.G", "SIE.G", 
                                                "comparison", "SRE.L vector", "SIE.L vector", 
                                                "SRE.G vector", "SIE.G vector", "CDE vector"))
  p3b$variable <- factor(p3b$variable, levels = c("baseline", 
                                                  "SRE.L", "SIE.L", "SRE.G", "SIE.G", 
                                                  "comparison", "SRE.L vector", "SIE.L vector", 
                                                  "SRE.G vector", "SIE.G vector", "CDE vector"))
  p4$variable <- factor(p4$variable, levels = c("baseline", 
                                                "SRE.L", "SIE.L", "SRE.G", "SIE.G", 
                                                "comparison", "SRE.L vector", "SIE.L vector", 
                                                "SRE.G vector", "SIE.G vector", "CDE vector"))
  p3b <- p3b[p3b$variable != "baseline", ]
  return(list(dataStandard,p2, p3b, p4))
}

named_group_split = function (.tbl, ...){
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), 
                                              sep = " / ")))
  grouped %>% group_split() %>% rlang::set_names(names)
}
####!!!-- END FUNCTIONS --!!!####
# grab the log by log diversity data
logDiversity = readRDS(file = "./Price-analysis/data/raw-data/logdiversity.rds") 
logDiversity_drop = readRDS(file = "./Price-analysis/data/raw-data/logdiversity_drop.rds")

# clean data
logList = logDiversity_drop %>% named_group_split(log_id) %>%
  purrr::map(~.x %>% 
               dplyr::select(log_id, xylo_Abundance, xylo_S,log10consumed, matches("^xylo\\d+")) %>%
               dplyr::mutate(massConsumed = (10^log10consumed),
                             perCapitaConsumed = massConsumed/xylo_Abundance) %>%
               dplyr::select(log_id, perCapitaConsumed, matches("^xylo\\d+")) %>%
               group_by(log_id, perCapitaConsumed) %>%
               pivot_longer(matches("^xylo\\d+"),names_to = 'species', values_to = 'abundance') %>%
               dplyr::mutate(abundance = replace_na(abundance, 0),
                             consumption = perCapitaConsumed*abundance) %>%
               ungroup)

# Determine the reference and treatment site order
pairwiseCombn = combn(names(logList), m = 2) %>% t() %>% 
  data.frame %>% 
  setNames(., nm = c('site1','site2')) %>%
  dplyr::mutate(id = 1:n()) %>%
  group_by(id) %>%
  pivot_longer(-id, names_to = 'siteNumber', values_to = 'log_id') %>%
  left_join(logDiversity_drop %>% dplyr::select(log_id, xylo_S), by = 'log_id') %>%
  left_join(logDiversity_drop %>% dplyr::select(log_id, log10consumed), by = 'log_id') %>% 
  left_join(logDiversity_drop %>% dplyr::select(log_id, xylo_Abundance), by = 'log_id') %>% 
  # internal helper function
  named_group_split(id) %>%
  # internal helper function
  lapply(., set_site_order) %>%
  bind_rows() %>% dplyr::select(id:log_id) %>%
  pivot_wider( names_from = siteNumber, values_from = log_id) %>%
  dplyr::select(site1, site2) %>%
  split(seq(nrow(.))) %>%
  lapply(., unlist)

#### Import data into priceTools format ----
# make_price_list -> internal helper function
priceDataSetup = purrr::map(pairwiseCombn, ~make_price_list(.x)) %>%
  purrr::map(~priceTools::data.setup(.x))  

#### Estimate the price partitions with `price.part()` ----
if(rerun){
pricePartitions = priceDataSetup %>%
  purrr::map(~priceTools::price.part(.x, sps.level = TRUE))
saveRDS(pricePartitions, "./Price-analysis/data/derived-data/pricePartitions.rds")
} else{
  pricePartitions = readRDS("./Price-analysis/data/derived-data/pricePartitions.rds")
}
#### Summarise the price components ----
# summarise_price -> internal helper function
priceEffectsDf = pricePartitions %>% purrr::map(~summarise_price(.x)) %>% 
  setNames(., lapply(pairwiseCombn, function(x) paste(unlist(x), collapse = "_"))) %>%
  bind_rows(.id = "compID") %>%
  remove_rownames() %>%
  dplyr::mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), NA,.x))) %>%
  tidyr::separate(col = compID, sep = "_", into = c('site1','site2'), remove = FALSE) %>%
  left_join(logDiversity_drop %>%
              dplyr::select(site1 = 'log_id', abundance_site1 = 'xylo_Abundance'), by = 'site1') %>%
  left_join(logDiversity_drop %>%
              dplyr::select(site2 = 'log_id', abundance_site2 = 'xylo_Abundance'), by = 'site2') %>%
  dplyr::mutate(deltaN = abundance_site1 - abundance_site2)

# summarise the price effects across pairwise comparisons
priceEffectsSumm = priceEffectsDf %>%
  dplyr::summarise(across(where(is.numeric), list(min = ~min(.x, na.rm = TRUE),
                                                  mean = ~mean(.x, na.rm = TRUE),
                                                  median = ~median(.x, na.rm = TRUE),
                                                  quant2.5 = ~quantile(.x, .025, na.rm = TRUE),
                                                  mad = ~mad(.x, na.rm = TRUE),
                                                  quant25 = ~quantile(.x, .25, na.rm = TRUE),
                                                  quant75 = ~quantile(.x, 0.75, na.rm = TRUE),
                                                  quant97.5 = ~quantile(.x, 0.975, na.rm = TRUE),
                                                  max = ~max(.x, na.rm = TRUE),
                                                  sd = ~sd(.x, na.rm = TRUE))))

#### Another way to partition out for plotting ----

priceEffectsDf2 = pricePartitions %>% purrr::map(~.x %>% pluck(1)) %>%
  setNames(., lapply(pairwiseCombn, function(x) paste(unlist(x), collapse = "_"))) %>%
  bind_rows(.id = 'compID')

# summarise these again #
priceEffectsSumm2 = priceEffectsDf2 %>%
  dplyr::summarise(across(where(is.numeric), list(min = ~min(.x, na.rm = TRUE),
                                                  mean = ~mean(.x, na.rm = TRUE),
                                                  median = ~median(.x, na.rm = TRUE),
                                                  quant2.5 = ~quantile(.x, .025, na.rm = TRUE),
                                                  mad = ~mad(.x, na.rm = TRUE),
                                                  quant25 = ~quantile(.x, .25, na.rm = TRUE),
                                                  quant75 = ~quantile(.x, 0.75, na.rm = TRUE),
                                                  quant97.5 = ~quantile(.x, 0.975, na.rm = TRUE),
                                                  max = ~max(.x, na.rm = TRUE),
                                                  sd = ~sd(.x, na.rm = TRUE))))

#### Arrange by difference in individual abundance
pairwiseCombnN = combn(names(logList), m = 2) %>% t() %>% 
  data.frame %>% 
  setNames(., nm = c('site1','site2')) %>%
  dplyr::mutate(id = 1:n()) %>%
  group_by(id) %>%
  pivot_longer(-id, names_to = 'siteNumber', values_to = 'log_id') %>%
  left_join(logDiversity_drop %>% dplyr::select(log_id, xylo_S, log10consumed, xylo_Abundance), by = 'log_id') %>%
  # left_join(logDiversity_drop %>% dplyr::select(log_id, log10consumed), by = 'log_id') %>% 
  # left_join(logDiversity_drop %>% dplyr::select(log_id, xylo_Abundance), by = 'log_id') %>% 
  # internal helper function
  named_group_split(id) %>%
  # internal helper function
  lapply(., set_site_order_n) %>%
  bind_rows() %>% dplyr::select(id:log_id) %>%
  pivot_wider( names_from = siteNumber, values_from = log_id) %>%
  dplyr::select(site1, site2) %>%
  split(seq(nrow(.))) %>%
  lapply(., unlist)

#### Import data into priceTools format ----
# make_price_list -> internal helper function
priceDataSetupN = purrr::map(pairwiseCombnN, ~make_price_list(.x)) %>%
  purrr::map(~priceTools::data.setup(.x))  

#### Estimate the price partitions with `price.part()` ----
if(rerun){
pricePartitionsN = priceDataSetupN %>%
  purrr::map(~priceTools::price.part(.x, sps.level = TRUE))
saveRDS(pricePartitionsN, "./Price-analysis/data/derived-data/pricePartitionsN.rds")
} else{
  pricePartitionsN = readRDS("./Price-analysis/data/derived-data/pricePartitionsN.rds")
}
# summarise_price -> internal helper function
priceEffectsNDf = pricePartitionsN %>% purrr::map(~summarise_price(.x)) %>% 
  setNames(., lapply(pairwiseCombnN, function(x) paste(unlist(x), collapse = "_"))) %>%
  bind_rows(.id = "compID") %>%
  remove_rownames() %>%
  dplyr::mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), NA,.x))) %>%
  tidyr::separate(col = compID, sep = "_", into = c('site1','site2'), remove = FALSE) %>%
  left_join(logDiversity_drop %>%
              dplyr::select(site1 = 'log_id', abundance_site1 = 'xylo_Abundance'), by = 'site1') %>%
  left_join(logDiversity_drop %>%
              dplyr::select(site2 = 'log_id', abundance_site2 = 'xylo_Abundance'), by = 'site2') %>%
  dplyr::mutate(deltaN = abundance_site1 - abundance_site2) %>%
  dplyr::mutate(across(c(SRE.L,SRE.G,SIE.L,SIE.G,CDE), ~(.x/deltaF)*100, .names = "{.col}rel"))

# extract just aggregated effects
priceEffectsNDf2 = pricePartitionsN %>% purrr::map(~.x %>% pluck(1)) %>%
  setNames(., lapply(pairwiseCombnN, function(x) paste(unlist(x), collapse = "_"))) %>%
  bind_rows(.id = 'compID')

# estimate the normalized effect for each component
priceEffectsNormalizedDf = pricePartitions %>% purrr::map(~.x %>% pluck(1)) %>%
  setNames(., lapply(pairwiseCombnN, function(x) paste(unlist(x), collapse = "_"))) %>%
  bind_rows(.id = 'compID') %>%
  dplyr::select(compID, SIE.L, SIE.G, SRE.L, SRE.G, CDE, x.func, y.func) %>%
  dplyr::mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), NA,.x))) %>%
  rowwise %>%
  dplyr::mutate(normalizer = max(abs(SRE.L), abs(SRE.G), abs(SIE.L), abs(SIE.G), abs(CDE))) %>%
  ungroup %>%
  dplyr::mutate(across(c(SIE.L, SIE.G, SRE.L, SRE.G, CDE), ~abs(.x)/normalizer)) %>%
  dplyr::select(-normalizer)

# create relative effect size for each price component
priceEffectsPropDf = priceEffectsDf2 %>%
  dplyr::mutate(across(c(SRE.L, SRE.G, SIE.L, SIE.G, CDE), ~(.x/x.func)*100))

# create df of aggregate effects: e.g., loss, gain, etc.
diffList = priceEffectsNDf2 %>%
  named_group_split(compID) %>%
  purrr::map(~process.data.price_cor(.x))

diffPriceDf = diffList %>% purrr::map(~.x %>% pluck(1)) %>%
  bind_rows(.id = 'compID') %>%
  left_join(priceEffectsNDf %>% dplyr::select(compID, deltaN)) %>%
  dplyr::mutate(SL = SIE.L+SRE.L,
                SG = SIE.G+SRE.G)

# create a component offset to force all effects to positive + 0.001
diffPriceSummDf = diffPriceDf %>%
  dplyr::select(compID, deltaN, SRE.L:CE) %>%
  dplyr::mutate(CDEmin = abs(min(CDE, na.rm = TRUE)),
                SRE.Lmin = abs(min(SRE.L, na.rm = TRUE)),
                SRE.Gmin = abs(min(SRE.G, na.rm = TRUE)),
                SIE.Lmin = abs(min(SIE.L, na.rm = TRUE)),
                SIE.Gmin = abs(min(SIE.G, na.rm = TRUE)),
                ) %>%
  dplyr::mutate(CDEtotal = CDE + (CDEmin+0.001),
                SRE.Ltotal = SRE.L + (SRE.Lmin+0.001),
                SRE.Gtotal = SRE.G + (SRE.Gmin + 0.001),
                SIE.Ltotal = SIE.L + (SIE.Lmin + 0.001),
                SIE.Gtotal = SIE.G + (SIE.Gmin + 0.001))

# relative (offset) price effects summary
relPriceDiffDf = diffPriceSummDf %>% dplyr::summarise(across(where(is.numeric),list(min = ~min(.x, na.rm = TRUE),
                                                                                    mean = ~mean(.x, na.rm = TRUE),
                                                                                    median = ~median(.x, na.rm = TRUE),
                                                                                    quant2.5 = ~quantile(.x, .025, na.rm = TRUE),
                                                                                    mad = ~mad(.x, na.rm = TRUE),
                                                                                    quant25 = ~quantile(.x, .25, na.rm = TRUE),
                                                                                    quant75 = ~quantile(.x, 0.75, na.rm = TRUE),
                                                                                    quant97.5 = ~quantile(.x, 0.975, na.rm = TRUE),
                                                                                    max = ~max(.x, na.rm = TRUE),
                                                                                    sd = ~sd(.x, na.rm = TRUE))))

# remove unneeded objects
rm(list = ls()[ls() %in% c("diffList","logDiversity","logDiversity_drop","logList","pairwiseCombn", "pairwiseCombnN","priceDataSetup","priceDataSetupN","priceEffectsDf2","priceEffectsDf","priceEffectsNDf","priceEffectsNDf2","priceEffectsSumm","priceEffectsSumm2","pricePartitions","pricePartitionsN","relPriceDiffDf")])