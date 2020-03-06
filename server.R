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
	v <- reactiveValues(file1=NULL, file2=NULL, file3=NULL, transcriptomics=NULL, trait=NULL, module=NULL, corr_type=NULL, df_output=data.frame())

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
	#	Gene-Module-Assignment information
	#	----------------------------------------------
	observeEvent(input$file3,{
		source("methods.R")
		v$file3=input$file3
		
		v$module=read.csv(input$file3$datapath,sep="\t",header=T,check.names=FALSE)

		my_mod=unique(v$module[,2])
		for(i in 1:length(my_mod)){
			nmb_genes=subset(v$module,v$module[,2]==my_mod[i])
			nmb_genes=dim(nmb_genes)[1]
			print(paste(my_mod[i],nmb_genes))
		}

		shinyalert("INFO", paste(dim(v$trait)[1],"gene-module assignments were uploaded!"), type = "info")
	})

	
	#	----------------------------------------------
	#	Calculation of transcriptomics and traits
	#	----------------------------------------------
	observeEvent(input$goButton,{
		source("methods.R")
		output$plot=renderPlot({
			par(oma=c(0,20,0,20))
			calc_cmp_transcriptomics_traits(v)
		})
		output$plot2=renderPlot({
		})
		output$plot3=renderPlot({
		})
	})
	
	#	----------------------------------------------
	#	Correlation of phenotypes
	#	----------------------------------------------	
	observeEvent(input$go_alpha2,{
		source("methods.R")
		output$plot=renderPlot({
			par(oma=c(0,20,0,20))
			L=make_corr(v,input$phen0,output)	
		})
		output$plot2=renderPlot({
		})
		output$plot3=renderPlot({
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
		output$plot2=renderPlot({
		})
		output$plot3=renderPlot({
		})
	})

	#	----------------------------------------------
	#	Regression analysis
	#	----------------------------------------------	
	observeEvent(input$go_alpha4,{
		source("methods.R")
		output$plot=renderPlot({
			regr_analysis(v,input$phen4,input$gene1)
		})
		output$plot2=renderPlot({
		})
		output$plot3=renderPlot({
		})
	})

	#	----------------------------------------------
	#	Linear model
	#	----------------------------------------------	
	observeEvent(input$go_alpha5,{
		source("methods.R")
		output$plot=renderPlot({
			res=linear_model(v,input$phen6,input$phen7,input$phen8,input$opt1,input$opt2,input$gene2,input$phen62,input$phen72,input$phen82)
			updateTextAreaInput(session, "phen5_area", label = "", value = res)

		})
		output$plot2=renderPlot({
		})
		output$plot3=renderPlot({
		})
	})

	#	----------------------------------------------
	#	Defining the correlation type
	#	----------------------------------------------
	observeEvent(input$corr_type,{
		v$corr_type=input$corr_type;
	})

	#	----------------------------------------------
	#	Multiple test correction
	#	----------------------------------------------
	observeEvent(input$multiple_test_correction,{
		v$multiple_test_correction=input$multiple_test_correction
	})

	#	----------------------------------------------
	#	Multiple test correction
	#	----------------------------------------------
	observeEvent(input$multiple_test_correction2,{
		v$multiple_test_correction2=input$multiple_test_correction2
	})
	
	observeEvent(input$significance_level,{
		v$significance_level=input$significance_level
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
	    }
	)
}
