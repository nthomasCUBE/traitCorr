# traitCorr (v0.8)

## TraitCorr - correlating gene expression with phenotypic data

The user interface TraitCorr allows to determine those genes that significantly correlate with a selected trait or among traits and provides various visualisation and analysis possibilities.


## Screenshots of the tool

correlation of transcriptome data with trait information

![traitCorr](https://github.com/nthomasCUBE/traitCorr/blob/master/pix/Figure1_V2.png)

## Literature

Please find our manuscript published in Gene Reports
[TraitCorr-Gene-Report-Nussbaumer et al. 2020](https://www.sciencedirect.com/science/article/pii/S2452014420300637)
and our preprint
[TraitCorr-Preprint-Nussbaumer et al., 2019](https://www.biorxiv.org/content/10.1101/557975v1)


## Installation

Following packages are needed, the Rscript "ReporteRs" has been replaced by R package "officer",
as suggestd by [davidgohel/ReporteRs](https://github.com/davidgohel/ReporteRs).
```
library(d3heatmap)
library(gplots)
library(officer)
library(VennDiagram)
library(scales)
library(shiny)
library(shinyalert)
library(shinyBS)
library(shinyjs)
library(shinythemes)
library(xlsx)
```

This is how to run it:

```
source("methods.R")
source("ui.R")
source("server.R")
shinyApp(ui,server)
```

To simplify the process of running traitcorr we added the method run_traitcorr that includes all dependencies.

Normally, I would copy the content of run_traitCorr.R into the R workbench and then
I first see whether R packages are still needed. If so, then the R shiny GUI should be displayed.
The input_data "traitCorr_input.zip" would have to get extracted and finally the data can be included
and afterwards possible to analyse the demo files.

```
run_traitcorr()
```



