# Emit the published mr.raps estimates for every exposure to outcome pair, and
# the loss function sensitivity on the primary analysis.
suppressWarnings(suppressMessages(library(mr.raps)))
H <- "results/stage7_matrix/harmonised"
ans <- c("biot->psa","tt->psa","shbg->psa","biot->prac","tt->prac","shbg->prac","biot->fg","psa->prac")
rows <- list()
for (an in ans) {
  m <- read.delim(file.path(H, paste0("harm_", gsub("->","_",an,fixed=TRUE), ".tsv")))
  r <- suppressWarnings(mr.raps(m$b_x, m$b_y, m$se_x, m$se_y,
                                over.dispersion = TRUE, loss.function = "huber"))
  rows[[length(rows)+1]] <- data.frame(analysis=an, beta=r$beta.hat, se=r$beta.se,
                                       tau2=r$tau2.hat, nsnp=nrow(m),
                                       package="mr.raps 0.2", settings="over.dispersion=TRUE, huber")
}
out <- do.call(rbind, rows)
out$lo <- out$beta - 1.96*out$se; out$hi <- out$beta + 1.96*out$se
out$p <- 2*pnorm(-abs(out$beta/out$se))
write.table(out, "results/stage8_crossvalidation/raps_reference.tsv", sep="\t", row.names=FALSE, quote=FALSE)

m <- read.delim(file.path(H,"harm_biot_psa.tsv"))
sens <- list()
for (lf in c("l2","huber","tukey")) for (od in c(TRUE, FALSE)) {
  r <- suppressWarnings(try(mr.raps(m$b_x,m$b_y,m$se_x,m$se_y,over.dispersion=od,loss.function=lf), silent=TRUE))
  if (!inherits(r,"try-error"))
    sens[[length(sens)+1]] <- data.frame(loss=lf, overdispersion=od, beta=r$beta.hat, se=r$beta.se,
      tau2=if(!is.null(r$tau2.hat)) r$tau2.hat else NA_real_)
}
s <- do.call(rbind, sens)
write.table(s, "results/stage8_crossvalidation/raps_loss_sensitivity.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("wrote reference RAPS for", nrow(out), "analyses and", nrow(s), "sensitivity rows\n")
print(out[,c("analysis","beta","se","tau2")], row.names=FALSE, digits=5)
