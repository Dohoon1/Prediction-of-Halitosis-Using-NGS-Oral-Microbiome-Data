# Salivary microbiome-based classification of concordant intra-oral halitosis

Code and data for the manuscript **“Salivary Microbiome-Based Classification of Concordant Intra-Oral Halitosis: *Porphyromonas gingivalis* as a Candidate Sentinel Taxon”** (Scientific Reports, revised submission).

Raw 16S rRNA gene (V3–V4) sequencing reads are available from the NCBI Sequence Read Archive under BioProject **PRJNA1148754**. The correspondence between SRA runs and study participants is given in `data/sample_metadata.csv` (and in Supplementary Table 1 of the manuscript).

Archived release: Zenodo, DOI 10.5281/zenodo.XXXXXXX (to be updated after the release is created).

## Repository structure

```
data/
  PCR_NGS_Data.xlsx       qPCR loads (16 species + total bacteria, copies/mL) and NGS ASV read counts for 104 participants
  sample_metadata.csv     participant ID, halitosis group, periodontal status, age, sex, smoking, oral hygiene,
                          mouth-air H2S / CH3SH / total VSCs (ppb), number of teeth/implants,
                          reads after QC, inclusion flag, SRA run and BioSample accessions
code/
  common.py               data loading and sample QC (shared by all scripts)
  01_analysis.py          all statistical analyses -> results/
  02_figures.py           all figures -> figures/
results/                  output tables (CSV/JSON) produced by 01_analysis.py
figures/                  Figures 1–7 and Supplementary Figures 1–3 (PNG and TIFF)
legacy_code/              scripts of the previous submission (single 80/20 split); kept for transparency
```

## How to run

```bash
pip install -r requirements.txt
cd code
python 01_analysis.py      # ~10 min on a laptop (permutation tests)
python 02_figures.py
```
All random procedures use fixed seeds; the outputs in `results/` and `figures/` are reproduced exactly.

## Sample inclusion

104 participants were enrolled (72 Non-Halitosis, 32 Halitosis). Sequencing failed for two participants (2S-16, 0 reads; 2S-20, 1 read); they are excluded by `common.py` (reads < 1,000), leaving **102 participants (72 Non-Halitosis, 30 Halitosis)**. These two participants have no SRA record.

## Correspondence between outputs and the manuscript

| Manuscript item | Output file(s) |
|---|---|
| Table 1 – age, sex, smoking, oral hygiene, VSC concentrations, periodontal status, number of teeth, *P. gingivalis* detection, total bacterial load, diversity | `results/R.json` (`T1`), `results/div.npy` |
| Table 2, Supplementary Table 2 – differential abundance (Wilcoxon + BH-FDR; sensitivity analyses) | `results/DA_full.csv` |
| Table 3 – cross-validated performance (repeated stratified 5-fold CV × 50, permutation tests) | `results/T3_cv.csv` |
| Table 4 – pairwise AUROC comparison (corrected resampled t-test, Holm) | `results/T4_pairs.csv` |
| Table 5 – demographic-augmented models (age, sex, smoking, oral hygiene) | `results/T5_demographic_models.csv` |
| Table 6 – linear regression for VSC concentrations (repeated 5-fold CV × 50) | `results/T6_vsc_regression.csv` |
| Supplementary Table 3 – qPCR comparisons and NGS–qPCR correlation | `results/qpcr.csv` |
| Supplementary Table 4 – *P. gingivalis* + random partner taxa | `results/null_pg_partners.csv` |
| Supplementary Table 5 – model coefficients (OR per 1-SD) | `results/FI.csv` |
| Supplementary Table 6 – models adjusted for periodontal status | `results/S_perio_models.csv` |
| Adjusted association of *P. gingivalis* with halitosis; within-stratum detection | `results/R.json` (`adj`, `adj_full`, `strat`) |
| Figures 4–5, Supplementary Figure 3 – within-group correlations | `results/corr_*.csv`, `results/corrp_*.csv` |
| Supplementary Figures 1–2 – composition, PCoA/PERMANOVA | `results/comp_*.csv`, `results/pcoa.npy`, `results/R.json` (`pcoa`) |

## Notes on naming

- Species-level assignments from partial 16S rRNA sequences are putative; ASV counts are summed by species label.
- The qPCR panel lists *Filifactor alocis* under its former name *Fusobacterium alocis*.
- Oral hygiene was rated good, moderate, or poor (no participant was rated poor).
- The prefixes 1S/2S/3S of participant IDs and the suffixes CON/GV/PO of SRA sample names denote the periodontal recruitment stratum (healthy/gingivitis/periodontitis), not the halitosis group.

## Software

Python 3.11; scikit-learn 1.8; statsmodels 0.15; SciPy 1.17; pandas; matplotlib.
