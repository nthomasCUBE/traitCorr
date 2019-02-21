library(d3heatmap)
library(shiny)
library(shinyalert)
library(shinyBS)
library(shinyjs)
library(shinythemes)

options(stringsAsFactors=FALSE)

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
		v$transcriptomics=read.csv(input$file1$datapath,sep="\t",header=T,check.names=FALSE)
		shinyalert("INFO", paste(dim(v$transcriptomics)[1],"transcripts from ",dim(v$transcriptomics)[2]," samples was uploaded!"), type = "info")
	})

	#	----------------------------------------------
	#	Trait information
	#	----------------------------------------------
	observeEvent(input$file2,{
		source("methods.R")
		v$file2=input$file2
		add_trait_information()
		v$trait=read.csv(input$file2$datapath,sep="\t",header=T,check.names=FALSE)
		shinyalert("INFO", paste(dim(v$trait)[1],"entries loaded for the trait information was uploaded!"), type = "info")
	})
	
	#	----------------------------------------------
	#	Calculation of transcriptomics and traits
	#	----------------------------------------------
	observeEvent(input$goButton,{
		source("methods.R")
		output$plot=renderPlot({
			calc_cmp_transcriptomics_traits(v)
		})
	})
	
	#	----------------------------------------------
	#	Correlation of phenotypes
	#	----------------------------------------------	
	observeEvent(input$go_alpha2,{
		source("methods.R")
		output$plot=renderPlot({
			print(paste0("phen0::",input$phen0))
			L=make_corr(v,input$phen0,output)	

		})
	})

	#	----------------------------------------------
	#	Comparison of different traits
	#	----------------------------------------------	
	observeEvent(input$go_alpha3,{
		source("methods.R")
		output$plot=renderPlot({
			cmp_traits(v,input$phen1,input$phen2,input$phen3)	
		})
	})
}
