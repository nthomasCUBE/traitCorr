#
# the following script is loading all libraries that are necessary
# to run traitCorr and finally open the R Shiny instance of traitCorr
#

packages_necessary=c("d3heatmap","gplots","officer","VennDiagram","scales","shiny","shinyalert","shinyBS","shinyjs","shinythemes","xlsx","plm")
for(x in 1:length(packages_necessary)){
  print(packages_necessary[x])
  library(packages_necessary[x],character.only=TRUE)
}

scripts_necessary=c("methods.R","ui.R","server.R")
for(x in 1:length(scripts_necessary)){
  print(scripts_necessary[x])
  source(scripts_necessary[x])
}

shinyApp(ui,server)
