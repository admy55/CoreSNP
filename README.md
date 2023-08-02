# Minimal coreSNP Selection
An Efficient Pipeline for Selecting Core Markers from Genome-wide SNP Datasets in Crops
<br>

We introduce the coreSNP pipeline, providing Shannon indices to evaluate the discrimination power of maker combinations, with the aim to select a minimal set of core SNPs essential for discrimination of large-scale sequenced samples. coreSNP is designed to use the compressed or uncompressed Variant Call Format (VCF) file as input to produce the core SNP sets results.

## Dependences
* Python >= 3.6
* Numpy
* [PLINK 1.9](https://www.cog-genomics.org/plink/)
* System platform: Linux, Windows, MacOS


## Usage
    select_coreSNP.py –-vcf <data.vcf>





![](images/pipeline.png)

<br>
