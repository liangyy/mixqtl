#!/bin/bash

#Rscript $MIXQTL/scripts/kk_run_mixqtl_conda.R \
Rscript $MIXQTL/scripts/mixqtl_conda.R \
-library_size $DATA/geuvadis.library_size.tsv.gz \
-variant_annotation $DATA/geuvadis_22_variant_annotation.txt.gz \
-window 1000000 \
-haplotype_1 $DATA/geuvadis_22_hap1.txt.gz \
-haplotype_2 $DATA/geuvadis_22_hap2.txt.gz \
-covariates $DATA/geuvadis.covariate.txt.gz \
-expression_annotation $DATA/gencode.v12.txt.gz \
-expression_total_count $DATA/expression_count.txt.gz \
-expression_count_1 $DATA/geuvadis.asc.h1.tsv.gz \
-expression_count_2 $DATA/geuvadis.asc.h2.tsv.gz \
-output mixqtl.txt.gz
