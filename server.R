library(d3heatmap)
library(officer)
library(shiny)
library(shinyalert)
library(shinyBS)
library(shinyjs)
library(shinythemes)

options(stringsAsFactors=FALSE)

server <- function(input, output, session)
{
	v <- reactiveValues(file1=NULL, file2=NULL, transcriptomics=NULL, trait=NULL, corr_type=NULL, df_output=data.frame())

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

	#	----------------------------------------------
	#	Defining the correlation type
	#	----------------------------------------------
	observeEvent(input$corr_type,{
		v$corr_type=input$corr_type;
	})

	output$download1 <- downloadHandler(
	    filename = function() {
	      paste("OTUs-traitCorr_report", ".csv", sep = "")
	    },
	    content = function(file) {
	      colnames(v$df_output)=c("correlation r","adj. p-value","Gene ID","trait")
	      write.csv2(v$df_output, file, row.names = FALSE)
	    }
	)

	output$download2 <- downloadHandler(
	    filename = function() {
	      paste("traitCorr_report", ".docx", sep = "")
	    },
	    content = function(file) {

	    my_doc=read_docx()
	    colnames(v$df_output)=c("correlation r","adj. p-value","Gene ID","trait")
	    body_add_table(my_doc,v$df_output, style = "table_template")
	    print(my_doc,file)
	    }
	)
}
