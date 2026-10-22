# Reproduce the authors' own R pipeline on the Dryad data, plus the unit-level
# alternatives the review needs. Writes tidy CSVs to handoff/r_results/.
D <- "handoff/dryad"; O <- "handoff/r_results"; dir.create(O, showWarnings = FALSE)
all      <- read.csv(file.path(D, "p1_all.csv"))
pairing  <- read.csv(file.path(D, "p1_pairings.csv"))
orders   <- read.csv(file.path(D, "p1_orders.csv"))
doglev   <- read.csv(file.path(D, "p1_dog_level.csv"))
longdat  <- read.csv(file.path(D, "p1_long_reconstructed.csv"))

res <- list(); k <- 0
addw <- function(id, label, x, y, data_sheet) {
  w <- suppressWarnings(wilcox.test(x, y, paired = TRUE))
  d <- x - y; d <- d[!is.na(d)]
  k <<- k + 1
  res[[k]] <<- data.frame(test_id = id, label = label, sheet = data_sheet, method = w$method,
                          statistic_name = "V", statistic = unname(w$statistic),
                          p_value = w$p.value, n_pairs = length(d),
                          n_nonzero = sum(d != 0), exact = grepl("exact", w$method, fixed = TRUE),
                          stringsAsFactors = FALSE)
}
# --- the three pooled tests, exactly as the authors wrote them (sheet 'all') ---
addw("W1", "habit1 vs habit3 (habituation)",     all$habit1, all$habit3, "all")
addw("W2", "habit4 vs novel (discrimination)",   all$habit4, all$novel,  "all")
addw("W3", "habit3 vs novel (dishabituation)",   all$habit3, all$novel,  "all")
addw("W4", "habit3 vs habit4 (format control)",  all$habit3, all$habit4, "all")
# --- by emotional pairing (sheet 'emotionalpairings', one row per dog) ---
addw("W5",  "JB habit4 vs novel", pairing$JBhabit4, pairing$JBnovel, "emotionalpairings")
addw("W6",  "JB habit3 vs novel", pairing$JBhabit3, pairing$JBnovel, "emotionalpairings")
addw("W11", "SB habit4 vs novel", pairing$SBhabit4, pairing$SBnovel, "emotionalpairings")
addw("W11b","SB habit3 vs novel", pairing$SBhabit3, pairing$SBnovel, "emotionalpairings")
addw("W12", "SJ habit4 vs novel", pairing$SJhabit4, pairing$SJnovel, "emotionalpairings")
addw("W14", "SJ habit3 vs novel", pairing$SJhabit3, pairing$SJnovel, "emotionalpairings")
# --- by presentation order (sheet 'specificorders', one row per dog) ---
addw("W7",  "J habituated / B novel: habit4 vs novel", orders$Jhabit4_vB, orders$Bnovel_vJ, "specificorders")
addw("W9",  "J habituated / B novel: habit3 vs novel", orders$Jhabit3_vB, orders$Bnovel_vJ, "specificorders")
addw("W8",  "B habituated / J novel: habit4 vs novel", orders$Bhabit4_vJ, orders$Jnovel_vB, "specificorders")
addw("W10", "B habituated / J novel: habit3 vs novel", orders$Bhabit3_vJ, orders$Jnovel_vB, "specificorders")
addw("O5",  "B habituated / S novel: habit4 vs novel", orders$Bhabit4_vS, orders$Snovel_vB, "specificorders")
addw("O6",  "B habituated / S novel: habit3 vs novel", orders$Bhabit3_vS, orders$Snovel_vB, "specificorders")
addw("O7",  "S habituated / B novel: habit4 vs novel", orders$Shabit4_vB, orders$Bnovel_vS, "specificorders")
addw("O8",  "S habituated / B novel: habit3 vs novel", orders$Shabit3_vB, orders$Bnovel_vS, "specificorders")
addw("W13", "S habituated / J novel: habit4 vs novel", orders$Shabit4_vJ, orders$Jnovel_vS, "specificorders")
addw("W15", "S habituated / J novel: habit3 vs novel", orders$Shabit3_vJ, orders$Jnovel_vS, "specificorders")
addw("O11", "J habituated / S novel: habit4 vs novel", orders$Jhabit4_vS, orders$Snovel_vJ, "specificorders")
addw("W16", "J habituated / S novel: habit3 vs novel", orders$Jhabit3_vS, orders$Snovel_vJ, "specificorders")
# --- the same pooled comparisons at dog level (mean of a dog's trials) ---
addw("W1_dog", "habit1 vs habit3, dog means",    doglev$habit1, doglev$habit3, "dog-level (43)")
addw("W2_dog", "habit4 vs novel, dog means",     doglev$habit4, doglev$novel,  "dog-level (43)")
addw("W3_dog", "habit3 vs novel, dog means",     doglev$habit3, doglev$novel,  "dog-level (43)")
addw("W4_dog", "habit3 vs habit4, dog means",    doglev$habit3, doglev$habit4, "dog-level (43)")
write.csv(do.call(rbind, res), file.path(O, "part1_wilcoxon.csv"), row.names = FALSE)

# --- paired t and sign test alongside, for the three-routes table -----------
rt <- list(); j <- 0
for (nm in list(c("W1","habit1","habit3"), c("W2","habit4","novel"), c("W3","habit3","novel"),
                c("W4","habit3","habit4"))) {
  for (lev in c("trial","dog")) {
    dd <- if (lev == "trial") all else doglev
    x <- dd[[nm[2]]]; y <- dd[[nm[3]]]; d <- x - y; d <- d[!is.na(d)]
    tt <- t.test(x, y, paired = TRUE)
    nz <- d[d != 0]
    st <- binom.test(sum(nz > 0), length(nz))
    j <- j + 1
    rt[[j]] <- data.frame(test_id = nm[1], level = lev, n = length(d), n_nonzero = length(nz),
                          mean_diff = mean(d), t_stat = unname(tt$statistic), t_df = unname(tt$parameter),
                          t_p = tt$p.value, ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
                          sign_pos = sum(nz > 0), sign_p = st$p.value)
  }
}
write.csv(do.call(rbind, rt), file.path(O, "part1_routes.csv"), row.names = FALSE)

# --- the six side-bias / season regressions --------------------------------
lmres <- list(); j <- 0
for (f in c("JBdiscrmdifference ~ JBside + Season", "JBdishabitdifference ~ JBside + Season",
            "SBdiscrmdifference ~ SBside + Season", "SBdishabitdifference ~ SBside + Season",
            "SJdiscrmdifference ~ SJside + Season", "SJdishabitdifference ~ SJside + Season")) {
  m <- lm(as.formula(f), data = pairing)
  s <- summary(m)
  co <- as.data.frame(s$coefficients)
  co$term <- rownames(co); co$formula <- f; co$n <- length(residuals(m))
  co$r_squared <- s$r.squared
  j <- j + 1; lmres[[j]] <- co
}
lmres <- do.call(rbind, lmres)
names(lmres)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
write.csv(lmres, file.path(O, "part1_sidebias_lm.csv"), row.names = FALSE)

# ===========================================================================
# Part 2
# ===========================================================================
sheets <- c("stressbehav","joybehav","owner","box","stranger","upright","sit","prone")
mf_sheets <- c("stressbehav","joybehav","owner","box","stranger")
pts <- paste0("X", seq(5, 60, 5), "s")

sw <- list(); kk <- list(); dn <- list(); dg <- list(); gl <- list()
for (i in seq_along(sheets)) {
  sh <- sheets[i]
  proc <- read.csv(file.path(D, paste0("p2_", sh, "_proc.csv")))
  rawd <- read.csv(file.path(D, paste0("p2_", sh, "_raw.csv")))
  s <- shapiro.test(proc$Frequency)
  sw[[i]] <- data.frame(sheet = sh, W = unname(s$statistic), p_value = s$p.value,
                        n = length(proc$Frequency),
                        route_taken = ifelse(sh == "box", "aov", "kruskal.test"))
  if (sh == "box") {
    a <- summary(aov(Frequency ~ Condition, data = proc))[[1]]
    kk[[i]] <- data.frame(sheet = sh, method = "aov", statistic_name = "F",
                          statistic = a[["F value"]][1], df1 = a[["Df"]][1], df2 = a[["Df"]][2],
                          p_value = a[["Pr(>F)"]][1], n = nrow(proc))
  } else {
    h <- kruskal.test(Frequency ~ Condition, data = proc)
    kk[[i]] <- data.frame(sheet = sh, method = "kruskal.test", statistic_name = "H",
                          statistic = unname(h$statistic), df1 = unname(h$parameter), df2 = NA,
                          p_value = h$p.value, n = nrow(proc))
  }
  # dog-level: each dog's mean and total across its sample points
  rawd$dog_mean  <- rowMeans(rawd[, pts], na.rm = TRUE)
  rawd$dog_total <- rowSums(rawd[, pts], na.rm = TRUE)
  hm <- kruskal.test(dog_mean ~ Odor, data = rawd)
  at <- summary(aov(dog_mean ~ Odor, data = rawd))[[1]]
  pm <- glm(dog_total ~ Odor, family = poisson, data = rawd)
  dg[[i]] <- data.frame(sheet = sh, n_dogs = nrow(rawd),
                        kw_H = unname(hm$statistic), kw_df = unname(hm$parameter), kw_p = hm$p.value,
                        aov_F = at[["F value"]][1], aov_p = at[["Pr(>F)"]][1],
                        pois_dev_p = anova(pm, test = "Chisq")[["Pr(>Chi)"]][2])
}
write.csv(do.call(rbind, sw), file.path(O, "part2_shapiro.csv"), row.names = FALSE)
write.csv(do.call(rbind, kk), file.path(O, "part2_omnibus.csv"), row.names = FALSE)
write.csv(do.call(rbind, dg), file.path(O, "part2_doglevel.csv"), row.names = FALSE)

# --- Dunn, both the package convention and the plain two-sided z -----------
dunn_manual <- function(x, g) {
  ok <- !is.na(x) & !is.na(g); x <- x[ok]; g <- factor(g[ok])
  N <- length(x); r <- rank(x)
  tb <- table(r); ties <- sum(tb^3 - tb)
  sig2 <- (N * (N + 1) / 12) - (ties / (12 * (N - 1)))
  lv <- levels(g); out <- list(); m <- 0
  for (a in 1:(length(lv) - 1)) for (b in (a + 1):length(lv)) {
    ia <- g == lv[a]; ib <- g == lv[b]
    na <- sum(ia); nb <- sum(ib)
    z <- (mean(r[ia]) - mean(r[ib])) / sqrt(sig2 * (1 / na + 1 / nb))
    m <- m + 1
    out[[m]] <- data.frame(comparison = paste(lv[a], "-", lv[b]), z = z, n1 = na, n2 = nb,
                           p_two_sided = 2 * pnorm(-abs(z)), p_half = pnorm(-abs(z)))
  }
  do.call(rbind, out)
}
dunn_sheets <- c("stressbehav","joybehav","owner","stranger","upright","prone")
dres <- list()
for (i in seq_along(dunn_sheets)) {
  sh <- dunn_sheets[i]
  proc <- read.csv(file.path(D, paste0("p2_", sh, "_proc.csv")))
  d <- dunn_manual(proc$Frequency, proc$Condition)
  d$sheet <- sh
  d$p_two_holm <- p.adjust(d$p_two_sided, method = "holm")
  d$p_half_holm <- p.adjust(d$p_half, method = "holm")
  d$p_two_bonf <- p.adjust(d$p_two_sided, method = "bonferroni")
  dres[[i]] <- d
}
write.csv(do.call(rbind, dres), file.path(O, "part2_dunn.csv"), row.names = FALSE)

# --- the five sex x condition GLMs, as written (no family = Gaussian) ------
for (i in seq_along(mf_sheets)) {
  sh <- mf_sheets[i]
  mf <- read.csv(file.path(D, paste0("p2_", sh, "_mf.csv")))
  m <- glm(Frequency_MF ~ Condition_MF * Sex_MF, data = mf)
  s <- summary(m)
  co <- as.data.frame(s$coefficients)
  co$term <- rownames(co); co$sheet <- sh; co$aic <- m$aic
  co$resid_df <- m$df.residual; co$n <- nrow(mf)
  co$dispersion <- s$dispersion
  gl[[i]] <- co
}
gl <- do.call(rbind, gl)
names(gl)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
write.csv(gl, file.path(O, "part2_glm_gaussian.csv"), row.names = FALSE)

# --- the same question at the dog level with a count family ----------------
gp <- list()
for (i in seq_along(mf_sheets)) {
  sh <- mf_sheets[i]
  rawd <- read.csv(file.path(D, paste0("p2_", sh, "_raw.csv")))
  rawd$dog_total <- rowSums(rawd[, pts], na.rm = TRUE)
  rawd$n_points <- rowSums(!is.na(rawd[, pts]))
  m <- glm(dog_total ~ Odor * Sex + offset(log(n_points)), family = quasipoisson, data = rawd)
  s <- summary(m)
  co <- as.data.frame(s$coefficients)
  co$term <- rownames(co); co$sheet <- sh; co$n_dogs <- nrow(rawd)
  co$dispersion <- s$dispersion; co$resid_df <- m$df.residual
  gp[[i]] <- co
}
gp <- do.call(rbind, gp)
names(gp)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
write.csv(gp, file.path(O, "part2_glm_dogcount.csv"), row.names = FALSE)

writeLines(paste("R", getRversion()), file.path(O, "r_version.txt"))
cat("done\n")
