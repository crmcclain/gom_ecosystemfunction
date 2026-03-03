suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(rlang)
  library(SIBER)
  library(MASS)      # mvrnorm
  library(Matrix)    # nearPD
  library(ggeffects) # ggpredict
})

# -----------------------------
# Choices
# -----------------------------
nreps <- 1999
set.seed(1)

pEll <- 0.95
min_n_total <- 10
min_n_group <- 5

# successful null reps (per log × metric) required for SES/p
min_null_success <- 50

# -----------------------------
# Helper: parse tree type from log ID (matches anywhere)
# -----------------------------
get_tree_from_log <- function(log_id) {
  if (grepl("L-G", log_id)) return("L-G")
  if (grepl("L-W", log_id)) return("L-W")
  NA_character_
}

# ============================================================
# SIBER helpers (community TNW + species SEAc + NI=SEAc/TNW)
# ============================================================

hullarea <- function(x, y) {
  ne <- length(x)
  abs(0.5 * ((x[1:(ne - 1)] %*% y[2:ne]) - (y[1:(ne - 1)] %*% x[2:ne])))
}

rini_from_kapow <- function(z) {
  tnw <- hullarea(z$bdry[[1]]$x, z$bdry[[1]]$y)
  
  seac <- NULL
  for (i in seq_along(z$owin.coords)) {
    temp.x <- z$owin.coords[[i]]$bdry[[1]]$x
    temp.y <- z$owin.coords[[i]]$bdry[[1]]$y
    seac <- rbind(seac, hullarea(temp.x, temp.y))
  }
  
  n.index <- NULL
  for (i in seq_along(z$owin.coords)) {
    temp.x <- z$owin.coords[[i]]$bdry[[1]]$x
    temp.y <- z$owin.coords[[i]]$bdry[[1]]$y
    n.index <- rbind(n.index, hullarea(temp.x, temp.y) / tnw)
  }
  
  id <- as.data.frame(unique(z$ell.coords$group))
  tnw <- as.data.frame(tnw)
  seac <- as.data.frame(seac)
  n.index <- as.data.frame(n.index)
  
  out <- cbind(id, tnw, seac, n.index)
  colnames(out) <- c("Species", "TNW", "SEAc", "Niche.Index")
  out
}

# ============================================================
# 1) Load isotope data and fit bivariate lognormal parameters
#    MVN on log(-d13C) and log(d15N)
# ============================================================

iso <- read.csv("/Users/Mark/Downloads/wdf_xylos.csv") %>%
  mutate(
    log = as.character(log),
    species = as.character(species),
    tree = vapply(log, get_tree_from_log, character(1)),
    d13C = as.numeric(d13C),
    d15N = as.numeric(d15N)
  ) %>%
  filter(
    !is.na(tree),
    is.finite(d13C), is.finite(d15N),
    d13C < 0,        # so -d13C > 0
    d15N > 0         # so log(d15N) is defined
  ) %>%
  mutate(
    Cpos = -d13C,
    Npos =  d15N,
    logC = log(Cpos),
    logN = log(Npos)
  )


fit_mvn2 <- function(logC, logN, ridge = 1e-6, min_n = 3) {
  x <- cbind(logC, logN)
  x <- x[complete.cases(x), , drop = FALSE]
  
  if (nrow(x) < min_n) {
    return(list(mu = c(NA_real_, NA_real_),
                Sigma = matrix(NA_real_, 2, 2),
                n = nrow(x)))
  }
  
  mu <- colMeans(x)
  Sigma <- stats::cov(x)
  
  # Guard: make Sigma positive definite
  ok_sigma <- all(is.finite(Sigma)) && isTRUE(all.equal(Sigma, t(Sigma)))
  if (!ok_sigma) {
    Sigma <- matrix(NA_real_, 2, 2)
  }
  
  if (any(!is.finite(Sigma)) || det(Sigma) <= 0) {
    Sigma <- as.matrix(Matrix::nearPD(Sigma)$mat)
  }
  Sigma <- Sigma + diag(ridge, 2)
  
  list(mu = mu, Sigma = Sigma, n = nrow(x))
}

# species × tree fit (preferred)
dist_sp_tree_bi <- iso %>%
  dplyr::group_by(species, tree) %>%
  dplyr::summarise(
    fit = list(fit_mvn2(logC, logN)),
    mu1 = fit[[1]]$mu[1],
    mu2 = fit[[1]]$mu[2],
    s11 = fit[[1]]$Sigma[1, 1],
    s12 = fit[[1]]$Sigma[1, 2],
    s21 = fit[[1]]$Sigma[2, 1],
    s22 = fit[[1]]$Sigma[2, 2],
    n_fit = fit[[1]]$n,
    .groups = "drop"
  ) %>%
  dplyr::select(-fit)

# species-only fallback
dist_sp_bi <- iso %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    fit = list(fit_mvn2(logC, logN)),
    mu1 = fit[[1]]$mu[1],
    mu2 = fit[[1]]$mu[2],
    s11 = fit[[1]]$Sigma[1, 1],
    s12 = fit[[1]]$Sigma[1, 2],
    s21 = fit[[1]]$Sigma[2, 1],
    s22 = fit[[1]]$Sigma[2, 2],
    n_fit = fit[[1]]$n,
    .groups = "drop"
  ) %>%
  dplyr::select(-fit)


get_bi_params_for_species_log <- function(species, log_id,
                                          dist_sp_tree_bi, dist_sp_bi,
                                          min_fit_n = 3) {
  
  tree <- get_tree_from_log(log_id)
  
  row <- dist_sp_tree_bi %>% filter(species == !!species, tree == !!tree)
  if (nrow(row) == 1 && is.finite(row$n_fit) && row$n_fit >= min_fit_n) {
    mu <- c(row$mu1, row$mu2)
    Sigma <- matrix(c(row$s11, row$s12,
                      row$s21, row$s22), 2, 2, byrow = TRUE)
    if (all(is.finite(mu)) && all(is.finite(Sigma))) return(list(mu = mu, Sigma = Sigma))
  }
  
  row2 <- dist_sp_bi %>% filter(species == !!species)
  if (nrow(row2) == 1 && is.finite(row2$n_fit) && row2$n_fit >= min_fit_n) {
    mu <- c(row2$mu1, row2$mu2)
    Sigma <- matrix(c(row2$s11, row2$s12,
                      row2$s21, row2$s22), 2, 2, byrow = TRUE)
    if (all(is.finite(mu)) && all(is.finite(Sigma))) return(list(mu = mu, Sigma = Sigma))
  }
  
  list(mu = c(NA_real_, NA_real_), Sigma = matrix(NA_real_, 2, 2))
}

# ---- simulate paired isotopes for one log from fitted MVN(logC, logN)
simulate_isotope_cloud_for_log <- function(count_vec_named, log_id,
                                           dist_sp_tree_bi, dist_sp_bi) {
  
  spp_present <- names(count_vec_named)[count_vec_named > 0]
  if (length(spp_present) == 0) {
    return(tibble(species = character(), d13C = numeric(), d15N = numeric(), log = character()))
  }
  ns <- as.integer(count_vec_named[spp_present])
  
  out <- purrr::map2_dfr(spp_present, ns, function(sp, n) {
    pars <- get_bi_params_for_species_log(sp, log_id, dist_sp_tree_bi, dist_sp_bi)
    if (!all(is.finite(pars$mu)) || !all(is.finite(pars$Sigma))) return(NULL)
    
    Z <- MASS::mvrnorm(n = n, mu = pars$mu, Sigma = pars$Sigma)
    
    Cpos <- exp(Z[, 1])
    Npos <- exp(Z[, 2])
    
    tibble(
      species = sp,
      d13C = -Cpos,
      d15N =  Npos,
      log  = log_id
    )
  })
  
  if (nrow(out) == 0) {
    out <- tibble(species = character(), d13C = numeric(), d15N = numeric(), log = character())
  }
  
  out
}

# ============================================================
# 2) Load abundance matrix and build assembly null
# ============================================================

#md2 <- readRDS("md.rds")
md2 <- as.matrix(md2)
md2 <- md2[order(rownames(md2)), , drop = FALSE]

logs_all <- colnames(md2)

meta_abund <- rowSums(md2)
meta_p <- meta_abund / sum(meta_abund)

S_obs <- colSums(md2 > 0)
N_obs <- colSums(md2)

# -----------------------------
# Choice: predictor diversity metric
# -----------------------------
pred_metric <- "Richness"  # c(Richness, Hill1, Hill2, N, Nmax, Pmax, S_eligible, rare_frac)

# -----------------------------
# Helpers: Hill numbers per log from abundance matrix
# -----------------------------
hill_numbers <- function(counts) {
  counts <- as.numeric(counts)
  N <- sum(counts, na.rm = TRUE)
  if (!is.finite(N) || N <= 0) {
    return(c(Richness = NA_real_, Hill1 = NA_real_, Hill2 = NA_real_))
  }
  
  p <- counts / N
  p <- p[p > 0 & is.finite(p)]
  
  Richness <- length(p)
  # Hill1 = exp(Shannon)
  Hill1 <- exp(-sum(p * log(p)))
  # Hill2 = inverse Simpson
  Hill2 <- 1 / sum(p^2)
  
  c(Richness = Richness, Hill1 = Hill1, Hill2 = Hill2)
}

extra_predictors <- function(counts, k = 5) {
  counts <- as.numeric(counts)
  N <- sum(counts, na.rm = TRUE)
  p <- if (N > 0) counts / N else rep(NA_real_, length(counts))
  p <- p[is.finite(p) & p > 0]
  
  S <- length(p)
  H <- if (S > 0) -sum(p * log(p)) else NA_real_
  Hill1 <- if (is.finite(H)) exp(H) else NA_real_
  Hill2 <- if (S > 0) 1 / sum(p^2) else NA_real_
  
  Nmax <- if (N > 0) max(counts, na.rm = TRUE) else NA_real_
  pmax <- if (S > 0) max(p, na.rm = TRUE) else NA_real_
  
  S_k <- sum(counts >= k, na.rm = TRUE)           # ellipse-eligible species count proxy
  rare_frac <- if (S > 0) mean(counts[counts > 0] < k) else NA_real_
  
  c(Richness = S, Hill1 = Hill1, Hill2 = Hill2,
    N = N, Nmax = Nmax, pmax = pmax,
    S_k = S_k, rare_frac = rare_frac)
}

k_group <- min_n_group  # use your same threshold

# allow for extra predictors
extra_predictors <- function(counts, k = 5) {
  counts <- as.numeric(counts)
  N <- sum(counts, na.rm = TRUE)
  p <- if (N > 0) counts / N else rep(NA_real_, length(counts))
  p <- p[is.finite(p) & p > 0]
  
  S <- length(p)
  H <- if (S > 0) -sum(p * log(p)) else NA_real_
  Hill1 <- if (is.finite(H)) exp(H) else NA_real_
  Hill2 <- if (S > 0) 1 / sum(p^2) else NA_real_
  
  Nmax <- if (N > 0) max(counts, na.rm = TRUE) else NA_real_
  pmax <- if (S > 0) max(p, na.rm = TRUE) else NA_real_
  
  S_k <- sum(counts >= k, na.rm = TRUE)           # ellipse-eligible species count proxy
  rare_frac <- if (S > 0) mean(counts[counts > 0] < k) else NA_real_
  
  c(Richness = S, Hill1 = Hill1, Hill2 = Hill2,
    N = N, Nmax = Nmax, pmax = pmax,
    S_k = S_k, rare_frac = rare_frac)
}

# make div_df
div_df <- tibble::tibble(log = logs_all) %>%
  mutate(stats = purrr::map(log, ~ extra_predictors(md2[, .x], k = k_group))) %>%
  mutate(
    Richness   = purrr::map_dbl(stats, ~ .x["Richness"]),
    Hill1      = purrr::map_dbl(stats, ~ .x["Hill1"]),
    Hill2      = purrr::map_dbl(stats, ~ .x["Hill2"]),
    N          = purrr::map_dbl(stats, ~ .x["N"]),
    Nmax       = purrr::map_dbl(stats, ~ .x["Nmax"]),
    pmax       = purrr::map_dbl(stats, ~ .x["pmax"]),
    S_eligible = purrr::map_dbl(stats, ~ .x["S_k"]),
    rare_frac  = purrr::map_dbl(stats, ~ .x["rare_frac"])
  ) %>%
  dplyr::select(-stats)

# Safety check (outdated)
#stopifnot(pred_metric %in% c("Richness", "Hill1", "Hill2"))

# ============================================================
# 3) SIBER metrics per log and species (SpLog) + log summaries
# ============================================================
 calc_siber_metrics_for_one_log <- function(points_log_df, pEll = 0.95,
                                           min_n_total = 10, min_n_group = 5) {
  
  # ---- NULL / type guard
  if (is.null(points_log_df) || !is.data.frame(points_log_df)) {
    return(tibble::tibble(
      TNW = NA_real_, mean_SEAc = NA_real_, mean_NicheIndex = NA_real_,
      n_groups_used = 0L, n_total = 0L
    ))
  }
  
  req <- c("species", "d13C", "d15N")
  if (!all(req %in% names(points_log_df))) {
    return(tibble::tibble(
      TNW = NA_real_, mean_SEAc = NA_real_, mean_NicheIndex = NA_real_,
      n_groups_used = 0L, n_total = nrow(points_log_df)
    ))
  }
  
  if (nrow(points_log_df) < min_n_total) {
    return(tibble::tibble(
      TNW = NA_real_, mean_SEAc = NA_real_, mean_NicheIndex = NA_real_,
      n_groups_used = 0L, n_total = nrow(points_log_df)
    ))
  }
  
  keep_groups <- points_log_df %>%
    dplyr::count(.data$species, name = "n") %>%
    dplyr::filter(.data$n >= min_n_group) %>%
    dplyr::pull(.data$species)
  
  df <- points_log_df %>%
    dplyr::filter(.data$species %in% keep_groups) %>%
    dplyr::filter(is.finite(.data$d13C), is.finite(.data$d15N))
  
  if (nrow(df) < min_n_total || dplyr::n_distinct(df$species) < 2) {
    return(tibble::tibble(
      TNW = NA_real_, mean_SEAc = NA_real_, mean_NicheIndex = NA_real_,
      n_groups_used = dplyr::n_distinct(df$species), n_total = nrow(df)
    ))
  }
  
  Y <- tibble::tibble(
    iso1  = df$d13C,
    iso2  = df$d15N,
    group = factor(df$species)
  )
  
  kap <- tryCatch(
    SIBER::siberKapow(Y, isoNames = c("iso1", "iso2"), group = "group", pEll = pEll),
    error = function(e) NULL
  )
  if (is.null(kap)) {
    return(tibble::tibble(
      TNW = NA_real_, mean_SEAc = NA_real_, mean_NicheIndex = NA_real_,
      n_groups_used = dplyr::n_distinct(df$species), n_total = nrow(df)
    ))
  }
  
  tab <- tryCatch(rini_from_kapow(kap), error = function(e) NULL)
  if (is.null(tab)) {
    return(tibble::tibble(
      TNW = NA_real_, mean_SEAc = NA_real_, mean_NicheIndex = NA_real_,
      n_groups_used = dplyr::n_distinct(df$species), n_total = nrow(df)
    ))
  }
  
  TNW <- suppressWarnings(as.numeric(tab$TNW[1]))
  mean_SEAc <- mean(as.numeric(tab$SEAc), na.rm = TRUE)
  mean_NI <- mean(as.numeric(tab$Niche.Index), na.rm = TRUE)
  
  tibble::tibble(
    TNW = TNW,
    mean_SEAc = mean_SEAc,
    mean_NicheIndex = mean_NI,
    n_groups_used = nrow(tab),
    n_total = nrow(df)
  )
}

calc_metrics_all_logs <- function(points_df, logs_all, pEll, min_n_total, min_n_group) {
  
  # split() drops empty logs; using [[L]] can return NULL => handled by function above
  spl <- split(points_df, points_df$log)
  
  out <- purrr::map_dfr(logs_all, function(L) {
    m <- calc_siber_metrics_for_one_log(
      points_log_df = spl[[L]],
      pEll = pEll, min_n_total = min_n_total, min_n_group = min_n_group
    )
    dplyr::mutate(m, log = L)
  })
  
  # put log first without relying on select() (avoids masking issues)
  out <- out[, c("log", setdiff(names(out), "log")), drop = FALSE]
  out
}

extract_species_metrics_for_one_log <- function(points_log_df, log_id, pEll = 0.95,
                                                min_n_total = 10, min_n_group = 5) {
  
  # NULL / type guard (fixes obs_points_spl[[L]] being NULL)
  if (is.null(points_log_df) || !is.data.frame(points_log_df)) return(NULL)
  
  req <- c("species", "d13C", "d15N")
  if (!all(req %in% names(points_log_df))) return(NULL)
  if (nrow(points_log_df) < min_n_total) return(NULL)
  
  keep_groups <- points_log_df %>%
    dplyr::count(.data$species, name = "n") %>%
    dplyr::filter(.data$n >= min_n_group) %>%
    dplyr::pull(.data$species)
  
  df <- points_log_df %>%
    dplyr::filter(.data$species %in% keep_groups) %>%
    dplyr::filter(is.finite(.data$d13C), is.finite(.data$d15N))
  
  if (dplyr::n_distinct(df$species) < 2) return(NULL)
  
  Y <- tibble::tibble(
    iso1  = df$d13C,
    iso2  = df$d15N,
    group = factor(df$species)
  )
  
  kap <- tryCatch(
    SIBER::siberKapow(Y, isoNames = c("iso1","iso2"), group = "group", pEll = pEll),
    error = function(e) NULL
  )
  if (is.null(kap)) return(NULL)
  
  tab <- tryCatch(rini_from_kapow(kap), error = function(e) NULL)
  if (is.null(tab)) return(NULL)
  
  tibble::tibble(
    log      = log_id,
    Species  = as.character(tab$Species),
    TNW      = as.numeric(tab$TNW),
    SEASpLog = as.numeric(tab$SEAc),          # <-- renamed from SEAc
    RNISpLog = as.numeric(tab$Niche.Index)
  )
}

# ============================================================
# 4) OBSERVED metrics (SIBER-based) from simulated "observed" communities
# ============================================================

obs_points <- map_dfr(logs_all, function(L) {
  counts <- setNames(as.integer(md2[, L]), rownames(md2))
  simulate_isotope_cloud_for_log(counts, L, dist_sp_tree_bi, dist_sp_bi)
})

obs_metrics <- calc_metrics_all_logs(obs_points, logs_all, pEll, min_n_total, min_n_group)

obs_TNW <- setNames(obs_metrics$TNW, obs_metrics$log)
obs_meanSEAc <- setNames(obs_metrics$mean_SEAc, obs_metrics$log)
obs_meanNI <- setNames(obs_metrics$mean_NicheIndex, obs_metrics$log)

# Species×log observed table (SpLog)
obs_points_spl <- split(obs_points, obs_points$log)
obs_spLog <- purrr::map_dfr(logs_all, function(L) {
  extract_species_metrics_for_one_log(
    points_log_df = obs_points_spl[[L]],
    log_id = L,
    pEll = pEll,
    min_n_total = min_n_total,
    min_n_group = min_n_group
  )
})
if (is.null(obs_spLog) || nrow(obs_spLog) == 0) {
  stop("obs_spLog is empty: species-level ellipses are failing everywhere.")
}

# ============================================================
# 5) NULL metrics (repeat nreps)
# ============================================================

null_TNW <- matrix(NA_real_, nrow = nreps, ncol = length(logs_all),
                   dimnames = list(paste0("rep", seq_len(nreps)), logs_all))
null_meanSEAc <- null_TNW
null_meanNI <- null_TNW

null_slopes_RNISpLog <- rep(NA_real_, nreps)  # RNISpLog ~ Richness
null_slopes_SEASpLog <- rep(NA_real_, nreps)  # SEAc (SpLog) ~ Richness

# helper: slope on species×log table joined to predictor (Richness/Hill1/Hill2)
slope_spLog_vs_predictor <- function(spLog_df, y, predictor = pred_metric) {
  if (is.null(spLog_df) || nrow(spLog_df) < 5) return(NA_real_)
  if (!all(c("log", y) %in% names(spLog_df))) return(NA_real_)
  
  dd <- spLog_df %>%
    dplyr::left_join(div_df %>% dplyr::select(log, all_of(predictor)), by = "log") %>%
    dplyr::rename(PRED = all_of(predictor)) %>%
    dplyr::filter(is.finite(.data[[y]]), is.finite(.data$PRED))
  
  if (nrow(dd) < 5) return(NA_real_)
  
  unname(coef(lm(reformulate("PRED", y), data = dd))["PRED"])
}


for (r in seq_len(nreps)) {
  
  # assemble null abundance matrix
  nullM <- matrix(0L, nrow(md2), ncol(md2), dimnames = dimnames(md2))
  for (j in seq_along(logs_all)) {
    nullM[, j] <- simulate_null_abund_one_log(S_obs[j], N_obs[j], meta_p, meta_abund)
  }
  
  # simulate isotope points under null
  null_points <- map_dfr(logs_all, function(L) {
    counts <- setNames(as.integer(nullM[, L]), rownames(nullM))
    simulate_isotope_cloud_for_log(counts, L, dist_sp_tree_bi, dist_sp_bi)
  })
  
  # log-level metrics
  null_metrics_r <- calc_metrics_all_logs(null_points, logs_all, pEll, min_n_total, min_n_group)
  null_TNW[r, ] <- null_metrics_r$TNW
  null_meanSEAc[r, ] <- null_metrics_r$mean_SEAc
  null_meanNI[r, ] <- null_metrics_r$mean_NicheIndex
  
  # species×log table for this rep -> slopes
  null_points_spl <- split(null_points, null_points$log)
  null_spLog_r <- purrr::map_dfr(logs_all, function(L) {
    extract_species_metrics_for_one_log(
      points_log_df = null_points_spl[[L]],
      log_id = L,
      pEll = pEll,
      min_n_total = min_n_total,
      min_n_group = min_n_group
    )
  })
  
  null_slopes_RNISpLog[r] <- slope_spLog_vs_predictor(null_spLog_r, "RNISpLog", predictor = pred_metric)
  null_slopes_SEASpLog[r] <- slope_spLog_vs_predictor(null_spLog_r, "SEASpLog", predictor = pred_metric)
  
  if (r %% 25 == 0) message("Done null rep ", r, "/", nreps)
}

# ============================================================
# 6) Summaries: SES / p-values (log-level metrics)
# ============================================================

summarize_vs_null <- function(obs_vec_named, null_mat, logs_all, min_null_success = 50) {
  
  obs_vec <- obs_vec_named[logs_all]
  
  n_finite <- apply(null_mat, 2, function(v) sum(is.finite(v)))
  null_mean <- apply(null_mat, 2, function(v) mean(v, na.rm = TRUE))
  null_sd   <- apply(null_mat, 2, function(v) sd(v, na.rm = TRUE))
  
  ok <- (n_finite >= min_null_success) &
    is.finite(null_sd) & (null_sd > 0) &
    is.finite(obs_vec)
  
  SES <- rep(NA_real_, length(logs_all))
  SES[ok] <- (obs_vec[ok] - null_mean[ok]) / null_sd[ok]
  
  z <- rep(NA_real_, length(logs_all))
  z[ok] <- SES[ok]  # same thing, just a familiar label
  
  p_lower <- rep(NA_real_, length(logs_all))
  p_upper <- rep(NA_real_, length(logs_all))
  p_two   <- rep(NA_real_, length(logs_all))
  
  for (j in seq_along(logs_all)) {
    if (!ok[j]) next
    v <- null_mat[, j]
    v <- v[is.finite(v)]
    p_lower[j] <- (sum(v <= obs_vec[j]) + 1) / (length(v) + 1)
    p_upper[j] <- (sum(v >= obs_vec[j]) + 1) / (length(v) + 1)
    p_two[j] <- min(1, 2 * min(p_lower[j], p_upper[j]))
  }
  
  tibble(
    log = logs_all,
    obs = as.numeric(obs_vec),
    null_mean = as.numeric(null_mean),
    null_sd = as.numeric(null_sd),
    n_null_finite = as.integer(n_finite),
    SES = as.numeric(SES),
    z = as.numeric(z),
    p_lower = as.numeric(p_lower),
    p_upper = as.numeric(p_upper),
    p_two = as.numeric(p_two),
    delta_obs_minus_exp = as.numeric(obs_vec - null_mean)
  )
}

res_TNW <- summarize_vs_null(obs_TNW, null_TNW, logs_all, min_null_success) %>%
  mutate(metric = "TNW (community 95% ellipse area)", axis = "2D (C,N)")
res_SEAc <- summarize_vs_null(obs_meanSEAc, null_meanSEAc, logs_all, min_null_success) %>%
  mutate(metric = "Mean SEAc (species 95% ellipse area)", axis = "2D (C,N)")
res_NI <- summarize_vs_null(obs_meanNI, null_meanNI, logs_all, min_null_success) %>%
  mutate(metric = "Mean Niche.Index (SEAc/TNW)", axis = "2D (C,N)")

res_siber <- bind_rows(res_TNW, res_SEAc, res_NI) %>%
  mutate(tree = vapply(log, get_tree_from_log, character(1)))

print(res_siber)

# ============================================================
# FINAL: three comparable null-slope histograms
#   1) TNW ~ Richness (log-level; uses res_siber SES table)
#   2) SEASpLog ~ Richness (species×log)
#   3) RNISpLog ~ Richness (species×log)
# Assumes you already ran:
#   - res_siber (with columns metric, SES, Richness) OR you can left_join Richness in
#   - null_TNW (nreps × logs matrix)
#   - null_slopes_SEASpLog, null_slopes_RNISpLog
#   - div_df (log, Richness)
# ============================================================

# ---- make sure res_siber has the chosen predictor column
if (!(pred_metric %in% names(res_siber))) {
  res_siber <- res_siber %>%
    dplyr::left_join(div_df %>% dplyr::select(log, all_of(pred_metric)), by = "log")
}


# -----------------------------
# Helper functions
# -----------------------------
get_obs_slope <- function(res_df, metric_label, response_col = "SES", predictor = pred_metric) {
  
  sub <- res_df %>%
    dplyr::filter(metric == metric_label,
                  is.finite(.data[[response_col]]),
                  is.finite(.data[[predictor]])) %>%
    dplyr::rename(PRED = all_of(predictor))
  
  if (nrow(sub) < 3) return(NA_real_)
  
  unname(coef(lm(reformulate("PRED", response_col), data = sub))["PRED"])
}

get_null_slopes <- function(null_mat, div_df, predictor = pred_metric, center_by = c("rep", "none")) {
  
  center_by <- match.arg(center_by)
  
  null_mat <- as.matrix(null_mat)
  stopifnot(!is.null(colnames(null_mat)))
  logs_all <- colnames(null_mat)
  nreps <- nrow(null_mat)
  
  # align predictor to log order in null_mat
  pred_vec <- div_df %>%
    dplyr::distinct(log, .data[[predictor]]) %>%
    dplyr::right_join(tibble::tibble(log = logs_all), by = "log") %>%
    dplyr::pull(.data[[predictor]])
  
  slopes <- rep(NA_real_, nreps)
  
  for (r in seq_len(nreps)) {
    
    df_r <- tibble::tibble(
      log = logs_all,
      null_value = as.numeric(null_mat[r, ]),
      PRED = pred_vec
    ) %>%
      dplyr::filter(is.finite(null_value), is.finite(PRED))
    
    if (nrow(df_r) < 3) next
    
    if (center_by == "rep") {
      df_r <- df_r %>% dplyr::mutate(null_centered = null_value - mean(null_value, na.rm = TRUE))
      y <- "null_centered"
    } else {
      y <- "null_value"
    }
    
    slopes[r] <- unname(coef(lm(reformulate("PRED", y), data = df_r))["PRED"])
  }
  
  slopes
}


plot_null_slope <- function(null_slopes, obs_slope, title, breaks = 30) {
  
  null_slopes <- null_slopes[is.finite(null_slopes)]
  
  if (!is.finite(obs_slope) || length(null_slopes) < 10) {
    stop("Need finite obs_slope and >=10 finite null slopes.")
  }
  
  ci <- quantile(null_slopes, probs = c(0.025, 0.975), na.rm = TRUE)
  
  max_abs <- max(abs(c(null_slopes, obs_slope)))
  pad <- 0.05 * max_abs
  if (!is.finite(pad) || pad == 0) pad <- 1
  xlim <- c(-(max_abs + pad), max_abs + pad)
  
  hist(null_slopes,
       breaks = breaks,
       xlim = xlim,
       col = "grey80",
       main = title,
       xlab = "Slope")
  
  abline(v = ci, lwd = 2, lty = 2, col = "black")  # 95% null interval
  abline(v = obs_slope, col = "red", lwd = 3)      # observed slope
}

# -----------------------------
# 1) TNW ~ Richness (log-level)
# -----------------------------
obs_slope_TNW <- get_obs_slope(
  res_siber,
  "TNW (community 95% ellipse area)",
  response_col = "SES",
  predictor = pred_metric
)

null_slopes_TNW <- get_null_slopes(null_TNW, div_df, predictor = pred_metric, center_by = "rep")

plot_null_slope(null_slopes_TNW, obs_slope_TNW,
                paste0("Total Niche Width (log) ~ ", pred_metric))


# -----------------------------
# 2) SEASpLog ~ Richness (species×log)
# -----------------------------
obs_slope_SEASpLog <- slope_spLog_vs_predictor(obs_spLog, "SEASpLog", predictor = pred_metric)

plot_null_slope(null_slopes_SEASpLog, obs_slope_SEASpLog,
                paste0("Standard Ellipse Area (species by log)", pred_metric))

# -----------------------------
# 3) RNISpLog ~ Richness (species×log)
# -----------------------------
obs_slope_RNISpLog <- slope_spLog_vs_predictor(obs_spLog, "RNISpLog", predictor = pred_metric)
plot_null_slope(null_slopes_RNISpLog, obs_slope_RNISpLog,
                paste0("Relative Niche Index (species by log)", pred_metric))

cor(div_df$Richness, div_df$Hill1, use = "complete.obs")
cor(div_df$Richness, div_df$Hill2, use = "complete.obs")
cor(div_df$Hill1, div_df$Hill2, use = "complete.obs")


