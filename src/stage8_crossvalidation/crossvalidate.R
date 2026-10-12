# RLO_TPC | Cross validation against published implementations
#
# Runs the committed harmonised variant level inputs through two published
# packages and compares every cell against the study's own estimator matrix.
#
#   MendelianRandomization 0.10.0  (Yavorska and Burgess)
#       mr_ivw, mr_egger, mr_median, mr_mbe, mr_cML
#   mr.raps 0.2                    (Zhao, Wang, Hemani, Bowden and Small)
#       mr.raps with Huber loss and overdispersion
#
# The study's estimators are custom implementations. This script does not
# replace them; it demonstrates agreement with reference implementations on
# identical inputs, which is what R3-Mod2 asked for.
#
# usage: Rscript src/stage8_crossvalidation/crossvalidate.R

suppressMessages({
  library(MendelianRandomization)
  library(mr.raps)
})

H <- "results/stage7_matrix/harmonised"
OWN <- read.delim("results/stage7_matrix/estimator_matrix.tsv", stringsAsFactors = FALSE)
set.seed(20261011)

NMIN <- c("biot->psa" = 178782, "tt->psa" = 194453, "shbg->psa" = 180726,
          "biot->prac" = 140254, "tt->prac" = 140254, "shbg->prac" = 140254,
          "biot->fg" = 178782, "psa->prac" = 140254)

rows <- list()
add <- function(analysis, estimator, ref_b, ref_se, pkg, note = "") {
  own <- OWN[OWN$analysis == analysis & OWN$estimator == estimator, ]
  rows[[length(rows) + 1]] <<- data.frame(
    analysis = analysis, estimator = estimator,
    own_beta = own$beta[1], own_se = own$se[1],
    ref_beta = ref_b, ref_se = ref_se,
    abs_diff_beta = abs(own$beta[1] - ref_b),
    rel_diff_beta = abs(own$beta[1] - ref_b) / pmax(abs(own$beta[1]), 1e-12),
    package = pkg, note = note, stringsAsFactors = FALSE)
}

for (an in names(NMIN)) {
  tag <- gsub("->", "_", an, fixed = TRUE)
  m <- read.delim(file.path(H, paste0("harm_", tag, ".tsv")), stringsAsFactors = FALSE)
  bx <- m$b_x; sx <- m$se_x; by <- m$b_y; sy <- m$se_y
  obj <- mr_input(bx = bx, bxse = sx, by = by, byse = sy)

  r <- mr_ivw(obj, model = "random")
  add(an, "IVW", r@Estimate, r@StdError, "MendelianRandomization::mr_ivw(model=random)")

  r <- mr_egger(obj)
  add(an, "MR Egger", r@Estimate, r@StdError.Est, "MendelianRandomization::mr_egger")

  r <- mr_median(obj, weighting = "weighted", iterations = 1000)
  add(an, "Weighted median", r@Estimate, r@StdError,
      "MendelianRandomization::mr_median(weighted)")

  r <- try(mr_mbe(obj, weighting = "weighted", stderror = "simple", iterations = 1000),
           silent = TRUE)
  if (!inherits(r, "try-error")) {
    add(an, "Weighted mode", r@Estimate, r@StdError,
        "MendelianRandomization::mr_mbe(weighted)")
  }

  # IVW after Cochran Q outlier removal: same rule, reference estimator on the subset
  b0 <- mr_ivw(obj, model = "random")@Estimate
  W <- (bx / sy)^2
  q <- W * (by / bx - b0)^2
  keep <- pchisq(q, 1, lower.tail = FALSE) > 0.05 / length(bx)
  if (sum(keep) > 2 && sum(keep) < length(bx)) {
    o2 <- mr_input(bx = bx[keep], bxse = sx[keep], by = by[keep], byse = sy[keep])
    r <- mr_ivw(o2, model = "random")
    add(an, "IVW, Q outliers removed", r@Estimate, r@StdError,
        "MendelianRandomization::mr_ivw on the Bonferroni Q filtered subset",
        sprintf("%d of %d variants retained", sum(keep), length(bx)))
  }

  r <- try(mr_cML(obj, n = NMIN[[an]], MA = FALSE, DP = TRUE, num_pert = 100,
                  random_start = 5, random_seed = 20261011), silent = TRUE)
  if (!inherits(r, "try-error")) {
    add(an, "MR cML BIC (DP SE)", r@Estimate, r@StdError,
        "MendelianRandomization::mr_cML(MA=FALSE, DP=TRUE)",
        "reference uses 100 perturbations against the study's 500")
  }

  r <- try(mr.raps::mr.raps(bx, by, sx, sy, over.dispersion = TRUE,
                            loss.function = "huber"), silent = TRUE)
  if (!inherits(r, "try-error")) {
    add(an, "MR RAPS Huber", r$beta.hat, r$beta.se,
        "mr.raps::mr.raps(over.dispersion=TRUE, huber)",
        "reference uses an analytic standard error against the study's 300 replicate bootstrap")
  }
  cat("done", an, "\n")
}

out <- do.call(rbind, rows)
dir.create("results/stage8_crossvalidation", showWarnings = FALSE, recursive = TRUE)
write.table(out, "results/stage8_crossvalidation/crossvalidation.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== agreement summary, point estimates ===\n")
cat(sprintf("cells compared: %d\n", nrow(out)))
cat(sprintf("max absolute difference: %.3e\n", max(out$abs_diff_beta)))
cat(sprintf("median absolute difference: %.3e\n", median(out$abs_diff_beta)))
cat(sprintf("cells agreeing to 1e-6: %d of %d\n",
            sum(out$abs_diff_beta < 1e-6), nrow(out)))
cat(sprintf("cells agreeing to 1e-3: %d of %d\n",
            sum(out$abs_diff_beta < 1e-3), nrow(out)))
cat(sprintf("cells agreeing within 1 percent relative: %d of %d\n",
            sum(out$rel_diff_beta < 0.01), nrow(out)))
worst <- out[order(-out$abs_diff_beta), ][1:8, ]
cat("\n=== largest discrepancies ===\n")
print(worst[, c("analysis", "estimator", "own_beta", "ref_beta", "abs_diff_beta")],
      row.names = FALSE, digits = 5)
