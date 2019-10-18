packages_necessary=c("d3heatmap","gplots","officer","VennDiagram","scales","shiny","shinyalert","shinyBS","shinyjs","shinythemes","xlsx")

source("methods.R")
source("ui.R")
source("server.R")

shinyApp(ui,server)
