# Plant- versus animal-based protein: an 8-week dietary substitution

Semi-targeted LC-MS/MS plasma metabolomics of the **SHAPES** cohort, from my PhD thesis
*Metabolomic phenotyping of different nutritional factors: focusing on meal timing and dietary
composition* (Department of Medicine, National University of Singapore, defended June 2026).

**Study design.** 8-week parallel-group randomised controlled trial, habitual protein-rich foods replaced with plant-based or animal-based products at matched labelled protein. The trial was designed and run by
collaborators, who provided plasma samples and de-identified clinical data for metabolomic
analysis.

This repository documents **methods and code**. It contains **no participant data**.

## Layout

```
analysis/    my analysis scripts for this cohort, plus shared_functions_SHAPES.R,
             the utility functions they call
data/        synthetic quant table documenting the expected input format
data-raw/    script that generates the synthetic table
```

## Assay and input

Agilent 1290 Infinity LC coupled to a 6495C triple quadrupole, SeQuant ZIC-cHILIC column,
dynamic MRM in positive and negative mode, roughly 300 to 500 analyte transitions spanning amino
acids, carnitines and acylcarnitines, nucleotides, organic acids and other polar metabolites.
Peak integration in [MRMkit](https://github.com/HuiSongLab/MRMkit).

MRMkit emits a tab-delimited quant table with three information columns, ten metric rows above
the sample rows (row 7 CoV, row 9 D-ratio, row 10 S/N), and pooled QC samples marked `BQC`.
`data/synthetic_quant_table.txt` is a synthetic table in that format, containing no real
measurements, so the reshaping logic can be checked without data.

## Quality control

| Metric | Threshold | Computed on |
|---|---|---|
| Signal-to-noise (S/N) | `>= 5` | uncorrected data |
| Coefficient of variation (CoV) | `<= 30%` | batch-corrected data |
| Dispersion ratio (D-ratio) | `<= 50%` | batch-corrected data |

An analyte is retained if `S/N >= 5 AND (CoV <= 30% OR D-ratio <= 50%)`. The CoV/D-ratio
criterion is deliberately a disjunction: CoV measures technical dispersion in absolute terms
while D-ratio expresses non-biological variance relative to biological variance, so an analyte
can be usable by one criterion and not the other and requiring both would discard signal.

Every surviving analyte is then inspected manually on four criteria, and analytes failing
manual review are excluded whether or not they passed the automated gate: peak shape, baseline
stability, retention-time consistency across samples, and agreement between subject samples and
pooled QC.

Intra-batch correction uses Gaussian kernel regression in MRMkit. Because the assay alternates
two HPLC columns to reduce cycle time, each column is treated as its own batch, so the batch
count passed to MRMkit is double the number of runs.

**Inter-batch grand-median correction.** In addition to MRMkit's intra-batch Gaussian kernel regression, this cohort uses an inter-batch correction implemented here: per analyte, the grand median across all batches is computed, then each batch is rescaled so its median matches, `corrected = raw * grand_median / batch_median`.

## Response measure

The response measure is a baseline-corrected metabolomic response (deltaMR) computed per analyte from week 0 and week 8 intensities. Because deltaMR can be negative and log2 ratios require positive values, an offset correction shifts all values into the positive range by adding a constant equal to the difference between 5 and the minimum negative deltaMR.

Throughout, the preprandial fasting metabolome is analysed separately from the postprandial or
post-intervention response. These answer different questions, and a difference present at
baseline is not evidence of a different response.

## Statistics


Mann-Whitney U on deltaMR between arms, with

```
log2FC = log2( median(adjusted deltaMR, arm A) / median(adjusted deltaMR, arm B) )
```

and median plus median absolute deviation for significant analytes. Benjamini-Hochberg FDR at adjusted p < 0.05, with analytes between 0.05 and 0.1 reported separately as trending rather than significant.

## Provenance and attribution

I analysed this cohort independently in R, and performed the LC-MS/MS assay for it myself on the laboratory's established platform. The utility functions my scripts call are in `analysis/shared_functions_SHAPES.R`. They **originate in earlier
analysis code from the laboratory**, inherited with the B&B project, and are not my own work. Where I
modified a function for this analysis its header is tagged "adapted from earlier lab code"; the rest
are as inherited.

See [ATTRIBUTION.md](ATTRIBUTION.md) for detail.

## Scope and limits

**No participant data is included.** Trial participant identifiers that appeared in the original
scripts have been redacted and replaced with `<..._REDACTED...>` placeholders, and absolute local
paths made relative. Affected files carry a header noting this.

**Results are not reported here.** Findings from this cohort are unpublished and belong to collaborative work in progress, so this repository describes how the analysis was done and not what it found.

**The scripts will not run end to end without the input data**, which is not distributed. The
shared laboratory functions they call are included as `analysis/shared_functions_SHAPES.R`; source that file
first, since in the original working environment these functions were loaded into the R session
directly rather than via `source()`.

## Dependencies

R, with `tidyverse`, `lme4`, `lmerTest`, `emmeans`, `qvalue` (Bioconductor), `mixOmics`,
`ggplot2`, `ggrepel`, `ggpubr`, `rstatix`, `car`, `gridExtra`. Peak integration and intra-batch
correction are performed externally in MRMkit; acquisition used Agilent MassHunter Workstation
10.0.127.

## Acknowledgements

The parent trial was **designed and conducted by Dr Darel Toh and colleagues**, who provided the
plasma samples and de-identified clinical data analysed here. The trial is registered with
ClinicalTrials.gov (NCT05446753) and was approved by the National Healthcare Group Domain Specific
Review Board, Singapore (reference 2022/00278).

I did not design or conduct the trial, recruit participants, or collect samples. My contribution
is the LC-MS/MS assay for this cohort, the metabolomic analysis, and the computational work in
this repository.

## A note on tooling

The documentation in this repository was drafted with AI assistance: this README,
[ATTRIBUTION.md](ATTRIBUTION.md), and the synthetic data generator in `data-raw/`.

The analysis code in `analysis/` was not. It is my own work, other than the laboratory's
pre-established utility functions identified in [ATTRIBUTION.md](ATTRIBUTION.md).

## Licence

Code I authored is released under the MIT Licence, see [LICENSE](LICENSE). `analysis/shared_functions_SHAPES.R` originates in earlier laboratory code and is not my own work.

## Contact

Dong Jiangwen. [Google Scholar](https://scholar.google.com/citations?user=NE3oTFMAAAAJ) ·
[LinkedIn](https://www.linkedin.com/in/jiangwen-dong)
