# Example Data

## Proteomics Data Source

The proteomics data in `example_proteome.csv` is from:

**Mulvey, J.F. et al.** (2026). An unbiased molecular characterisation of peripartum cardiomyopathy hearts identifies mast cell chymase as a new diagnostic candidate. *Molecular & Cellular Proteomics*, Volume 0, Issue 0, 101510.

## Sample Metadata

The sample metadata in `old_sample_meta.csv` was derived from the proteomics data patterns using the `generate_sample_metadata.R` script in the parent directory. This metadata is for demonstration purposes and is not from the original study.

The metadata includes:
- `sample_id` - Sample identifiers
- `group` - Control or heart failure
- `Sex` - Male or Female (randomly assigned)
- `Age` - Age in years (derived from PC2 + noise for demonstration)
- `BMI` - Body mass index (random)
- `Ejection fraction (%)` - Cardiac function measure (strongly associated with group and PC1)

The ejection fraction values are biologically realistic and correlated with the proteomics data structure to enable meaningful demonstration of PC-metadata association testing.
