# MrGWASSS
MrGWASSS cleans and visualizes GWAS summary statistics



## Dependencies

```txt
Depends: 
    BioThingsClient.R (>= 0.99.1),
    devtools,
    dplyr,
    ggplot2,
    ggpubr,
    knitr,
    latex2exp,
    purrr,
    writexl
Remotes:  
    biothings/BioThingsClient.R
```



## Usage

#### Load the package

```R
# Load the package
if(!require("MrGWASSS")){
  devtools::install_github("Broccolito/MrGWASSS")
  library(MrGWASSS)
}
```



#### Load GWAS summary statistics into the workspace

Given the load_data function arguments:

```R
load_data(
  file_name = "FHS_EA_MRS_5e8_snplist.txt",
  marker_name_column = "SNPID",
  chr_column = "CHR",
  pos_column = "POS",
  ref_column = "NEA",
  alt_column = "EA",
  pvalue_column = "p.value",
  delimiter = " "
)
```



To load the summary statistics from SAIGE, use:

```R
data = load_data("HAPO_AA_Ogtt_diff_chr_MERGED.txt", 
                 by_marker = TRUE, 
                 marker_name_column = "SNPID", 
                 pvalue_column = "p.value", 
                 delimiter = " ")
```



To load the summary statistics from METAL, use:

```R
data = load_data("METAANALYSIS1.TBL",
                 by_marker = TRUE,
                 marker_name_column = "MarkerName",
                 pvalue_column = "P.value",
                 delimiter = "\t")
```



#### Threshold the summary statistics

```R
data_1e5 = threshold_data(data, 1e-5)
```



#### Annotate the stats

```R
data_annotated = annotate_data(data_1e5)
```



#### Make Manhattan plots

```R
plot_manhattan(data)
```



#### Make Q-Q plot

```R
plot_qq(data)
```

