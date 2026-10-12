suppressMessages({library(cause); library(dplyr); library(readr)})
args <- commandArgs(trailingOnly=TRUE); f1<-args[1]; f2<-args[2]; tag<-args[3]
NSUB <- 100000   # nuisance-parameter subset, reduced for a 3 GB machine
set.seed(20261012)
X1 <- read_tsv(f1, col_types="cddcc"); X2 <- read_tsv(f2, col_types="cddcc")
X <- gwas_merge(X1, X2, snp_name_cols=c("snp","snp"), beta_hat_cols=c("beta","beta"),
                se_cols=c("se","se"), A1_cols=c("A1","A1"), A2_cols=c("A2","A2"))
rm(X1,X2); gc()
cat(tag,"merged:",nrow(X),"\n", file=stderr())
set.seed(20261012)
varlist <- sample(X$snp, size=min(NSUB,nrow(X)), replace=FALSE)
params <- est_cause_params(X, varlist)
cat(tag,"params done\n", file=stderr())
saveRDS(params, paste0(tag,"_params.rds"))
# LD pruning using the clumped instrument-style approach on p-value
X$p1 <- 2*pnorm(-abs(X$beta_hat_1/X$seb1))
cand <- X %>% filter(p1 < 1e-3)
write_tsv(data.frame(SNP=cand$snp, P=cand$p1), paste0(tag,"_cand.txt"))
cat(tag,"candidates:",nrow(cand),"\n", file=stderr())
