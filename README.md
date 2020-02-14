# traitCorr (v0.7)

## TraitCorr - correlating gene expression with phenotypic data

The user interface TraitCorr allows to determine those genes that significantly correlate with a selected trait or among traits and provides various visualisation and analysis possibilities.


## Screenshots of the tool

correlation of transcriptome data with trait information

![traitCorr](https://github.com/nthomasCUBE/traitCorr/blob/master/pix/Figure1_V2.png)

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

```
run_traitcorr()
```



