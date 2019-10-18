packages_necessary=c("d3heatmap","gplots","officer","VennDiagram","scales","shiny","shinyalert","shinyBS","shinyjs","shinythemes","xlsx")
for(x in 1:length(packages_necessary)){
  print(packages_necessary[x])
  library(packages_necessary[x],character.only=TRUE)
}

source("methods.R")
source("ui.R")
source("server.R")

shinyApp(ui,server)
