# This script runs GAM models for all Price components
# These can be rerun or loaded from the script `price_work.R`

# Run models of CDE patterns
CDE_n_brm = brms::brm(brms::bf(log10(CDEtotal)~s(deltaN)),
                      data = diffPriceSummDf,
                      family = gaussian(),
                      cores = 4, seed = 42,
                      iter = 7000, warmup = 5000, thin = 10,
                      control = list(adapt_delta = 0.995),
                      save_pars = save_pars(all = TRUE),
                      file = "./Price-analysis/data/derived-data/models/CDE_n_brm.rds",
                      file_refit = 'on_change')
# add in criterion
CDE_n_brm = add_criterion(CDE_n_brm, c('loo'))

# compare null model
CDE_null = brms::brm(brms::bf(log10(CDEtotal)~1),
                     data = diffPriceSummDf,
                     family = gaussian(),
                     cores = 4, seed = 42,
                     iter = 7000, warmup = 5000, thin = 10,
                     control = list(adapt_delta = 0.995),
                     save_pars = save_pars(all = TRUE),
                     file = "./Price-analysis/data/derived-data/models/CDE_null.rds",
                     file_refit = 'on_change')
# add criterion
CDE_null = add_criterion(CDE_null, c('loo'))

# summaries
summary(CDE_n_brm)
summary(CDE_null)
loo::loo_compare(CDE_n_brm, CDE_null)

#estimated change
CDE_n_ceff <- conditional_effects(CDE_n_brm, 'deltaN')
plot(CDE_n_ceff)

## SRE.LG
# Run models of SRE.G patterns
SREG_n_brm = brms::brm(brms::bf(log10(SRE.Gtotal)~s(deltaN)),
                       data = diffPriceSummDf,
                       family = gaussian(),
                       cores = 4, seed = 42,
                       iter = 7000, warmup = 5000, thin = 10,
                       control = list(adapt_delta = 0.995),
                       save_pars = save_pars(all = TRUE),
                       file = "./Price-analysis/data/derived-data/models/SRE.G_n_brm.rds",
                       file_refit = 'on_change')
SREG_n_brm = add_criterion(SREG_n_brm, c('loo'))

SREG_null = brms::brm(brms::bf(log10(SRE.Gtotal)~1),
                      data = diffPriceSummDf,
                      family = gaussian(),
                      cores = 4, seed = 42,
                      iter = 7000, warmup = 5000, thin = 10,
                      control = list(adapt_delta = 0.995),
                      save_pars = save_pars(all = TRUE),
                      file = "./Price-analysis/data/derived-data/models/SRE.G_null.rds",
                      file_refit = 'on_change')
SREG_null = add_criterion(SREG_null, c('loo'))

loo::loo_compare(SRE.G_n_brm, SRE.G_null)

SREG_n_ceff <- conditional_effects(SREG_n_brm, 'deltaN')
plot(SREG_n_ceff)

SREL_n_brm = brms::brm(brms::bf(log10(SRE.Ltotal)~s(deltaN)),
                       data = diffPriceSummDf,
                       family = gaussian(),
                       cores = 4, seed = 42,
                       iter = 7000, warmup = 5000, thin = 10,
                       control = list(adapt_delta = 0.995),
                       save_pars = save_pars(all = TRUE),
                       file = "./Price-analysis/data/derived-data/models/SRE.L_n_brm.rds",
                       file_refit = 'on_change')
SREL_n_brm = add_criterion(SREL_n_brm, c('loo'))

summary(SREL_n_brm)

SREL_null = brms::brm(brms::bf(log10(SRE.Ltotal)~1),
                      data = diffPriceSummDf,
                      family = gaussian(),
                      cores = 4, seed = 42,
                      iter = 7000, warmup = 5000, thin = 10,
                      control = list(adapt_delta = 0.995),
                      save_pars = save_pars(all = TRUE),
                      file = "./Price-analysis/data/derived-data/models/SRE.L_null.rds",
                      file_refit = 'on_change')
SREL_null = add_criterion(SREL_null, c('loo'))
loo::loo_compare(SRE.L_n_brm, SRE.L_null)

SREL_n_ceff <- conditional_effects(SREL_n_brm, 'deltaN')
plot(SREL_n_ceff)

diffPriceSummDf %>%
  dplyr::filter(SIE.Ltotal > 1 ) %>%
  ggplot()+
  geom_point(aes(x = deltaN, y = log10(SIE.Ltotal))) +
  geom_smooth(aes(x = deltaN, y = log10(SIE.Ltotal)))

SIEG_n_brm = brms::brm(brms::bf(log10(SIE.Gtotal)~s(deltaN, k = 6)),
                       data = diffPriceSummDf %>% dplyr::filter(SIE.Gtotal > 200),
                       family = gaussian(),
                       cores = 4, seed = 42,
                       iter = 9000, warmup = 7000, thin = 10,
                       control = list(adapt_delta = 0.995),
                       save_pars = save_pars(all = TRUE),
                       file = "./Price-analysis/data/derived-data/models/SIE.G_n_brm.rds",
                       file_refit = 'on_change')
SIEG_n_brm = add_criterion(SIEG_n_brm, c('loo'))

SIEG_null = brms::brm(brms::bf(log10(SIE.Gtotal)~1),
                      data = diffPriceSummDf %>% dplyr::filter(SIE.Gtotal > 200),
                      family = gaussian(),
                      cores = 4, seed = 42,
                      iter = 9000, warmup = 7000, thin = 10,
                      control = list(adapt_delta = 0.995),
                      save_pars = save_pars(all = TRUE),
                      file = "./Price-analysis/data/derived-data/models/SIE.G_null.rds",
                      file_refit = 'on_change')
SIEG_null = add_criterion(SIEG_null, c('loo'))

loo::loo_compare(SIE.G_n_brm, SIE.G_null)

SIEG_n_ceff <- conditional_effects(SIEG_n_brm, 'deltaN')
plot(SIEG_n_ceff)

SIEL_n_brm = brms::brm(brms::bf(log10(SIE.Ltotal)~s(deltaN)),
                       data = diffPriceSummDf %>% dplyr::filter(SIE.Ltotal > 1),
                       family = gaussian(),
                       cores = 4, seed = 42,
                       iter = 7000, warmup = 5000, thin = 10,
                       control = list(adapt_delta = 0.995),
                       save_pars = save_pars(all = TRUE),
                       file = "./Price-analysis/data/derived-data/models/SIE.L_n_brm.rds",
                       file_refit = 'on_change')
SIEL_n_brm = add_criterion(SIEL_n_brm, c('loo'))

SIEL_null = brms::brm(brms::bf(log10(SIE.Ltotal)~1),
                      data = diffPriceSummDf%>% dplyr::filter(SIE.Ltotal > 1),
                      family = gaussian(),
                      cores = 4, seed = 42,
                      iter = 7000, warmup = 5000, thin = 10,
                      control = list(adapt_delta = 0.995),
                      save_pars = save_pars(all = TRUE),
                      file = "./Price-analysis/data/derived-data/models/SIE.L_null.rds",
                      file_refit = 'on_change')
SIEL_null = add_criterion(SIEL_null, c('loo'))

loo::loo_compare(SIEL_n_brm, SIEL_null)

SIEL_n_ceff <- conditional_effects(SIEL_n_brm, 'deltaN')
plot(SIEL_n_ceff)