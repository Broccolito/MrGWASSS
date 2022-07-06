# MrGWASSS
MrGWASSS Cleans and Visualizes GWAS Summary Statistics



### Example

```R
### Example

# Install the package
library(devtools)
install_github("Broccolito/MrGWASSS")

# Load the package
library("MrGWASSS")

# Load GWAS summary statistics into the workspace
data = load_data("HAPO_AA_Ogtt_diff_chr_MERGED.txt")

# Threshold the stats
data_1e5 = threshold_data(data, 1e-5)

# Annotate the stats
data_annotated = annotate_data(data_1e5)

# Make Manhattan plot
plot_manhattan(data)

# Make Q-Q plot
plot_qq(data)
```

