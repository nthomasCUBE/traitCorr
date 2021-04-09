library(shiny)
library(shinyalert)
library(shinyBS)
library(shinyjs)
library(shinythemes)
library(officer)
library(xlsx)


options(stringsAsFactors=FALSE)
options(shiny.maxRequestSize = 50*1024^2)

ui <- fluidPage(  

tags$head(
	tags$style(HTML("
	.shiny-output-error {
	visibility: hidden;
}
body {
	#background-color: #23443333;
}
body, label, input, button, select { 
	font-family: 'Arial';
}"))
  ), 
  theme = shinytheme("united"),  useShinyjs(), useShinyalert(), 
	sidebarLayout(
		sidebarPanel(
		tabsetPanel(id = "tabset",
		tabPanel("TraitCorr - finding significantly correlating traits",
			fileInput("file1", "Transcriptome (expression data)", multiple = TRUE, accept = c("text/text", ".txt")),
			fileInput("file2", "Trait information", multiple = TRUE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
			fileInput("file3", "Gene-Module-Assignment", multiple = TRUE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
			actionButton("goButton", "Analyse dataset!")
		)
		)
		),
		mainPanel(
			useShinyjs(),
			plotOutput(outputId = "plot"),
			plotOutput(outputId = "plot2"),
			plotOutput(outputId = "plot3")
		)
	)
)