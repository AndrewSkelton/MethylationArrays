# Cartilage Methylation Processes
This repository is a snapshot of cleaned R Markdown scripts used in the analysis of Illumina HumanMethylation 450K bead chips for a cohort of cartilage samples. The primary disease of study is Osteoarthritis

## Repository Structure
├── 01_Preprocessing.Rmd
├── 02_Modelling_Differential_Methyl.Rmd
├── 03_cAge_Regression.Rmd
├── 04_Spline_Fits.Rmd
├── 05_Hovarth_Methylation_Clock.Rmd
├── Hovarth
│   ├── AdditionalFile3.csv
│   ├── CpG_Query.xlsx
│   ├── Horvath_MethyAge_Normalisation_Functions.R
│   ├── MethylationDataExample55.csv
│   ├── NORMALIZATION.R
│   ├── datMiniAnnotation.csv
│   └── probeAnnotation21kdatMethUsed.csv
├── README.md
└── Reference
    ├── HumanMethylation450_hg19_bowtie_multimap.txt
    └── crossreact-probes-Illumina450k.csv

* 01_Preprocessing.Rmd - Array preprocessing procedure
* 02_Modelling_Differential_Methyl.Rmd - Modelling differential methylation between groups
* 03_cAge_Regression.Rmd - Regression of chronological age (cAge), stratified by condition
* 04_Spline_Fits.Rmd - Spline fits of cAge with 2 degrees of freedom
* 05_Hovarth_Methylation_Clock.Rmd - Prediction of Methylation Age (mAge), per Hovarth's method.
* Hovarth - Reference datasets for Hovarth's method
* Reference - Reference datasets for probe filtering
