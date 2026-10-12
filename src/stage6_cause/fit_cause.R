suppressMessages({library(cause); library(dplyr); library(readr)})
args <- commandArgs(trailingOnly=TRUE); f1<-args[1]; f2<-args[2]; tag<-args[3]; prun<-args[4]
set.seed(20261012)
X1 <- read_tsv(f1, col_types="cddcc"); X2 <- read_tsv(f2, col_types="cddcc")
X <- gwas_merge(X1, X2, snp_name_cols=c("snp","snp"), beta_hat_cols=c("beta","beta"),
                se_cols=c("se","se"), A1_cols=c("A1","A1"), A2_cols=c("A2","A2"))
rm(X1,X2); gc()
params <- readRDS(paste0(tag,"_params.rds"))
top <- read.table(prun, header=TRUE)$SNP
top <- intersect(top, X$snp)
cat(tag," variants for CAUSE: ",length(top),"\n", file=stderr())
res <- cause(X=X, variants=top, param_ests=params, force=TRUE)
saveRDS(res, paste0(tag,"_cause.rds"))
s <- summary(res, ci_size=0.95)
cat("=== ",tag," ===\n")
print(s$tab)
print(s$p)
q <- res$elpd
print(q)
