# mixqtl

We provide a unified framework that combines total and allele-specific (at the gene level) read counts and is scalable to large studies with thousands of samples, such as GTEx. 

We develop computationally efficient methods that integrate allele-specific and total read counts to improve cisQTL discovery, fine-mapping, and genetic prediction of the transcriptome. Our log-linear approximation allows scaling-up of computations well-beyond a few hundreds of samples, providing an alternative to current total read count-based methods.

This repository hosts the R implementation of mixQTL, mixFine, and mixPred.


# How to use

We analyze one gene at a time. 
For a given gene, the following input data is required

* Total read count
* Allele-specific read count for haplotype 1 and 2 respectively
* Genotype of cis-window for haplotype 1 and 2 respectively
* Library size
* Covariates (optional)

# Tutorial for applying mixQTL
Find detailed instructions on how to run mixQTL [here](https://github.com/hakyimlab/mixqtl/wiki/Example-and-tutorial)
