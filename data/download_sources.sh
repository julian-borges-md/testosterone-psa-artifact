#!/usr/bin/env bash
# Retrieves every public input and verifies SHA256 against docs/DATA_CHECKSUMS.txt.
# Raw files land in work/ (git ignored). Requires curl and sha256sum.
set -euo pipefail
mkdir -p work/p2 work/tools && cd work
B=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
get(){ curl -sSL --retry 5 -C - -o "$1" "$2"; }
get biot_m.tsv.gz        $B/GCST90012001-GCST90013000/GCST90012103/GCST90012103_buildGRCh37.tsv.gz
get shbg_m.tsv.gz        $B/GCST90012001-GCST90013000/GCST90012109/GCST90012109_buildGRCh37.tsv.gz
get tt_m.tsv.gz          $B/GCST90012001-GCST90013000/GCST90012113/GCST90012113_buildGRCh37.tsv.gz
get psa.tsv.gz           $B/GCST90461001-GCST90462000/GCST90461907/GCST90461907.tsv.gz
get pca_practical.tsv.gz $B/GCST006001-GCST007000/GCST006085/harmonised/29892016-GCST006085-EFO_0001663-build37.f.tsv.gz
get pca_finngen.gz       https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_PROSTATE_EXALLC.gz
get p2/bmi.gz            $B/GCST002001-GCST003000/GCST002783/All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz
get p2/igf1.h.tsv.gz     $B/GCST90019001-GCST90020000/GCST90019511/harmonised/GCST90019511.h.tsv.gz
get p2/fi.tsv.gz         $B/GCST90002001-GCST90003000/GCST90002238/GCST90002238_buildGRCh37.tsv.gz
sha256sum -c ../docs/DATA_CHECKSUMS.txt
echo "LD reference: download http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz and extract EUR.* into work/tools/"
echo "PLINK 1.90b7.2: https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip into work/tools/"
