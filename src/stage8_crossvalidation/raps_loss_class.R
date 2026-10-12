# MR RAPS loss function sensitivity for the CLASS of statements that rest on the
# MR RAPS row, not only for the primary analysis.
#
# Four statements in the main text rest on the same MR RAPS row: the
# bioavailable testosterone count, the SHBG count, the contrast count and the
# stated p upper bound. A loss sensitivity tabulated for the primary analysis
# alone leaves the other three undetermined, so this script computes all three
# PSA exposures under each loss and the bioavailable minus total contrast
# within each loss. The alternative counts are then measured rather than
# described as unknown.
suppressWarnings(suppressMessages(library(mr.raps)))
H <- "results/stage7_matrix/harmonised"
EXP <- c("biot->psa", "tt->psa", "shbg->psa")
LOSS <- c("l2", "huber", "tukey")

fit <- function(an, lf) {
  m <- read.delim(file.path(H, paste0("harm_", gsub("->", "_", an, fixed = TRUE), ".tsv")))
  r <- suppressWarnings(try(mr.raps(m$b_x, m$b_y, m$se_x, m$se_y,
                                    over.dispersion = TRUE, loss.function = lf), silent = TRUE))
  if (inherits(r, "try-error")) return(NULL)
  data.frame(analysis = an, loss = lf, beta = r$beta.hat, se = r$beta.se,
             tau2 = if (!is.null(r$tau2.hat)) r$tau2.hat else NA_real_, nsnp = nrow(m))
}

rows <- list()
for (an in EXP) for (lf in LOSS) {
  f <- fit(an, lf)
  if (!is.null(f)) rows[[length(rows) + 1]] <- f
}
E <- do.call(rbind, rows)
E$lo <- E$beta - 1.96 * E$se
E$hi <- E$beta + 1.96 * E$se
E$p <- 2 * pnorm(-abs(E$beta / E$se))
E$excludes_null <- (E$lo > 0 & E$hi > 0) | (E$lo < 0 & E$hi < 0)
write.table(E, "results/stage8_crossvalidation/raps_loss_by_exposure.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# contrast within each loss, same construction as the committed contrast table:
# difference of the two point estimates with the standard errors combined in
# quadrature, the two GWAS being independent samples
con <- list()
for (lf in LOSS) {
  b <- E[E$analysis == "biot->psa" & E$loss == lf, ]
  t <- E[E$analysis == "tt->psa" & E$loss == lf, ]
  if (nrow(b) != 1 || nrow(t) != 1) next
  d <- b$beta - t$beta
  sd <- sqrt(b$se^2 + t$se^2)
  con[[length(con) + 1]] <- data.frame(loss = lf, biot = b$beta, tt = t$beta,
                                       difference = d, se = sd, z = d / sd,
                                       p = 2 * pnorm(-abs(d / sd)),
                                       both_negative = (b$beta < 0 & t$beta < 0))
}
C <- do.call(rbind, con)
write.table(C, "results/stage8_crossvalidation/raps_loss_contrast.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n== MR RAPS by exposure and loss ==\n")
print(E[, c("analysis", "loss", "beta", "se", "lo", "hi", "p", "excludes_null")],
      row.names = FALSE, digits = 4)
cat("\n== bioavailable minus total contrast, within loss ==\n")
print(C[, c("loss", "difference", "se", "z", "p")], row.names = FALSE, digits = 4)
