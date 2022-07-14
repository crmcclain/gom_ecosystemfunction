# load packages, import and clean data
source("./Price-analysis/code/01_price-data-load.R")

priceEffectsNormalizedDfSumm = priceEffectsNormalizedDf %>%
  ungroup %>%
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

priceEffectsNormalizedDfSumm %>%
  dplyr::select(matches("mean")) %>% data.frame

priceEffectsNormalizedDfSumm %>%
  dplyr::select(matches("sd")) %>% data.frame

priceEffectsPropDfSumm = priceEffectsPropDf %>% ungroup %>%
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

priceEffectsPropDfSumm %>%
  dplyr::select(matches("mean")) %>% data.frame
priceEffectsPropDfSumm %>%
  dplyr::select(matches("median")) %>% data.frame
priceEffectsPropDfSumm %>%
  dplyr::select(matches("mad")) %>% data.frame

## boxplots of price effects ----
# normalized price

priceNorm_plot = ggplotGrob(priceEffectsNormalizedDf %>%
  dplyr::select(compID, SRE.L, SRE.G, SIE.L, SIE.G, CDE) %>%
  pivot_longer(-compID, names_to = "component", values_to = 'value') %>%
  dplyr::mutate(component = factor(component, levels = c("SRE.L", "SRE.G", "SIE.L", "SIE.G", "CDE"))) %>%
  ggplot()+
  # geom_hline(aes(yintercept = 0), linetype = 'dashed', size = 1.1)+
  geom_violin(aes(x = component, y = value, fill = component), width = 0.7, alpha = 0.5, adjust = 0.5 )+
  scale_y_continuous(name = "Relative Importance", limits = c(0,1), expand = c(0.01,0.01))+
  scale_x_discrete()+
  annotate('text', label = 'A', x = 0, y = 1, hjust = 0)+
  viridis::scale_fill_viridis(discrete = TRUE, begin = 1,end =0)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        legend.position = 'none'));dev.off();grid.draw(priceNorm_plot)

priceProp_plot = ggplotGrob(priceEffectsPropDf %>%
                              dplyr::select(compID, SRE.L, SRE.G , SIE.L, SIE.G , CDE) %>%
                              pivot_longer(-compID, names_to = "component", values_to = 'value') %>%
                              dplyr::mutate(component = factor(component, levels = c("SRE.L", "SRE.G", "SIE.L", "SIE.G", "CDE"))) %>%
                              ggplot()+
                              geom_hline(aes(yintercept = 0), linetype = 'dashed', size = 1.1)+
                              geom_violin(aes(x = component, y = value, fill = component), width = 0.7, alpha = 0.5, adjust = 0.7)+
                              scale_y_continuous(name = "Contribution to mass loss (%)", expand = c(0.01,0.01))+
                              coord_cartesian(ylim = c(-200,200))+
                              scale_x_discrete()+
                              annotate('text', label = 'B', x = 0, y = Inf, hjust = 0, vjust = 1)+
                              viridis::scale_fill_viridis(discrete = TRUE, begin = 1,end =0)+
                              theme_minimal()+
                              theme(axis.title.x = element_blank(),
                                    legend.position = 'none'));dev.off();grid.draw(priceProp_plot)

# GAM models ----
# code to run GAMs in brms is `run-GAM-models.R`
# skip rerunning and import all model objects previously saved
load_models = TRUE
if(load_models){
  files = list.files("./Price-analysis/data/derived-data/models/","brm|null.rds", full.names = TRUE)
  model_names = gsub(".rds", "",sapply(strsplit(files, "/"),"[", 6))
  model_list = lapply(files, readRDS) %>% setNames(., nm = model_names)
  # load model objects to Global
  list2env(model_list, globalenv())
}

# model effects and comparisons -----
# conditional effects objects
SREL_n_ceff <- conditional_effects(SRE.L_n_brm, 'deltaN')
SREG_n_ceff <- conditional_effects(SRE.G_n_brm, 'deltaN')
SIEL_n_ceff <- conditional_effects(SIE.L_n_brm, 'deltaN')
SIEG_n_ceff <- conditional_effects(SIE.G_n_brm, 'deltaN')
CDE_n_ceff <- conditional_effects(CDE_n_brm, 'deltaN')
# null model comparisonsl
loo::loo_compare(SRE.L_n_brm, SRE.L_null)
loo::loo_compare(SRE.G_n_brm, SRE.G_null)
loo::loo_compare(SIE.L_n_brm, SIE.L_null)
loo::loo_compare(SIE.G_n_brm, SIE.G_null)
loo::loo_compare(CDE_n_brm, CDE_null)
# individual model criterion
loo::loo(SRE.L_n_brm)
loo::loo(SRE.G_n_brm)
loo::loo(SIE.L_n_brm)
loo::loo(SIE.G_n_brm)
loo::loo(CDE_n_brm)

# prediction plots
CDE_ceff_plot = 
  data.frame(
    variable = "CDE",
    deltaN = CDE_n_ceff$deltaN$deltaN,
    estimate = CDE_n_ceff$deltaN$estimate__,
    upper = CDE_n_ceff$deltaN$upper__,
    lower = CDE_n_ceff$deltaN$lower__
    )

SRE.L_ceff_plot = 
  data.frame(
    variable = "SRE.L",
    deltaN = SREL_n_ceff$deltaN$deltaN,
    estimate = SREL_n_ceff$deltaN$estimate__,
    upper = SREL_n_ceff$deltaN$upper__,
    lower = SREL_n_ceff$deltaN$lower__
  )
SRE.G_ceff_plot = 
  data.frame(
    variable = "SRE.G",
    deltaN = SREG_n_ceff$deltaN$deltaN,
    estimate = SREG_n_ceff$deltaN$estimate__,
    upper = SREG_n_ceff$deltaN$upper__,
    lower = SREG_n_ceff$deltaN$lower__
  )
SIE.L_ceff_plot = 
  data.frame(
    variable = "SIE.L",deltaN = SIEL_n_ceff$deltaN$deltaN,
    estimate = SIEL_n_ceff$deltaN$estimate__,
    upper = SIEL_n_ceff$deltaN$upper__,
    lower = SIEL_n_ceff$deltaN$lower__
  )
SIE.G_ceff_plot = 
  data.frame(
    variable = "SIE.G",
    deltaN = SIEG_n_ceff$deltaN$deltaN,
    estimate = SIEG_n_ceff$deltaN$estimate__,
    upper = SIEG_n_ceff$deltaN$upper__,
    lower = SIEG_n_ceff$deltaN$lower__
  )

# Back correction of all price effects to natural scales
priceCeffDf = bind_rows(SRE.L_ceff_plot,
                      SRE.G_ceff_plot,
                      SIE.L_ceff_plot,
                      SIE.G_ceff_plot,
                      CDE_ceff_plot
                      ) %>%
  dplyr::mutate(variable = factor(variable, levels = c("SRE.L", "SRE.G", "SIE.L", "SIE.G", "CDE"))) %>%
  left_join(diffPriceSummDf %>% dplyr::select(matches('min')) %>%
              pivot_longer(everything(),names_to = "variable", values_to = "correction") %>%
              distinct %>%
              dplyr::mutate(variable = gsub("min", "", variable),
                            variable = factor(variable, levels = c("SRE.L", "SRE.G", "SIE.L", "SIE.G", "CDE")))) %>%
  dplyr::mutate(across(c(estimate, upper, lower), ~(10^(.x))-(correction+0.001)))
  

PriceN_plot = priceCeffDf %>%
  ggplot()+
  geom_ribbon(aes(x = (deltaN), ymin = (lower), ymax = (upper), fill = variable), alpha = 0.4)+
  geom_line(aes(x = (deltaN), y = (estimate), color = variable, linetype = variable), size = 1.1)+
  scale_y_continuous(name = expression("Effect Contribution to mass loss (%)"),
                     expand = c(0.001,0.001))+
  scale_x_continuous(name = expression(""*Delta~N), limits = c(0,NA),expand = c(0.01,0.01))+
  coord_cartesian(ylim = c(-300, 500))+
  annotate('text', label = 'C', x = 0, y = Inf, hjust = 0, vjust = 1)+
  scale_linetype_manual(values = c("dotted", "solid", "dotted", rep("solid", 2)))+
  viridis::scale_fill_viridis(discrete = TRUE, begin = 1, end= 0)+
  viridis::scale_color_viridis(discrete = TRUE, begin = 1, end = 0)+
  theme(legend.title = element_blank(),
        legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.direction = 'horizontal');PriceN_plot

#
pdf("./Price-analysis/figures/Price_nEffects.pdf",width = 7, height = 7)
grid.arrange(grobs = list(priceNorm_plot, priceProp_plot ,PriceN_plot),
                        layout_matrix = rbind(cbind(1,2),
                                              cbind(3,3)), ncol = 2)
dev.off()

#### End of code ####