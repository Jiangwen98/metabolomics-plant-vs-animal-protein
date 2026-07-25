# Attribution and provenance

Laboratory analysis code is collaborative and cumulative. Utilities get written, inherited and
rewritten across projects, and projects themselves get handed over. This file records that
lineage plainly.

## This cohort

I analysed this cohort independently in R, and performed the LC-MS/MS assay for it myself on the laboratory's established platform. The plotting functions my scripts call build on earlier analysis code from the laboratory, which I adapted for this analysis. They are in `analysis/shared_functions_SHAPES.R`.

## Analysis functions

`analysis/shared_functions_SHAPES.R` holds the utility functions my scripts call: quality-control filtering,
transformation, scaling, outlier handling, volcano and time-course plotting.

The utility functions my scripts call are in `analysis/shared_functions_SHAPES.R`. These are the **Drum
Laboratory's pre-established analysis code**, which predates my project and came to me with the B&B
handover. They are **not my own work**.

I have deliberately not attributed them to an individual. They were developed within the laboratory,
with contributions and guidance from more than one member over time, so naming a single author would
be a claim I cannot verify. Credit belongs to the laboratory and to the people who built them.

Where I modified a function for the requirements of this analysis, its header is tagged "adapted from
earlier lab code"; the rest are as inherited.

The scripts do not `source()` this file, because in the original working environment the functions
were loaded into the R session directly. Source it first to run anything.

## My own contribution

The thesis attribution statement reads:

> "The plasma samples and original study designs of the three cohorts were provided by our
> collaborators. The LC-MS/MS assay panel and analytical method used throughout this thesis are an
> established platform of the Drum Laboratory. Within that platform, the thesis author
> independently performed the LC-MS/MS assay of the SHAPES cohort and carried out the processing of
> all raw metabolomics data (chromatographic peak integration and quality control), the downstream
> statistical analysis, and the interpretation of the results."

Consistent with that, my contribution here is the analysis pipeline itself: the choice and
sequencing of methods, the quality-control strategy and its thresholds, the response measure, the
statistical models and contrasts, and the interpretation. Calling shared utility functions no more
makes the analysis someone else's than calling a package does.

Some regions of the scripts also carry comments from other lab members, reflecting discussion and
review during the work.

## Acknowledgements

The parent trial was **designed and conducted by Dr Darel Toh and colleagues**, who provided the
plasma samples and de-identified clinical data analysed here. The trial is registered with
ClinicalTrials.gov (NCT05446753) and was approved by the National Healthcare Group Domain Specific
Review Board, Singapore (reference 2022/00278).

I did not design or conduct the trial, recruit participants, or collect samples. My contribution
is the LC-MS/MS assay for this cohort, the metabolomic analysis, and the computational work in
this repository.

## Third-party software

- **MRMkit**, peak integration and Gaussian-kernel-regression batch correction. Third-party
  software, not authored here.
- **Agilent MassHunter Workstation 10.0.127**, data acquisition.
- R packages as listed in the README, each under its own licence.
