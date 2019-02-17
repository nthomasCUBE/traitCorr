library(d3heatmap)
library(shiny)
library(shinyalert)
library(shinyBS)
library(shinyjs)
library(shinythemes)

server <- function(input, output, session)
{
	v <- reactiveValues(file1=NULL, file2=NULL, transcriptomics=NULL, trait=NULL)

	#	----------------------------------------------
	#	Transcriptomics
	#	----------------------------------------------
	observeEvent(input$file1,{
		source("methods.R")
		v$file1=input$file1
		add_transcriptomics()
		v$transcriptomics=read.csv(input$file1$datapath,sep="\t",header=T)
		shinyalert("INFO", "transcriptome was uploaded!", type = "info")
	})

	#	----------------------------------------------
	#	Trait information
	#	----------------------------------------------
	observeEvent(input$file2,{
		source("methods.R")
		v$file2=input$file2
		add_trait_information()
		v$trait=read.csv(input$file2$datapath,sep="\t",header=T)
		shinyalert("INFO", "trait information was uploaded!", type = "info")
	})
}
