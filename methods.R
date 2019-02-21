library(d3heatmap)
library(gplots)
library(VennDiagram)
library(scales)
library(xlsx)

options(stringsAsFactors=FALSE)

add_transcriptomics=function(){
	print(c("INFO|add_transcriptomcis"))
}

add_trait_information=function(){
	print(c("INFO|add_trait_information"))
}

calc_cmp_transcriptomics_traits=function(v){
	print(c("INFO|transcriptomics VS traits"))
	cn=c("---",colnames(v$trait))
	appendTab(inputId = "tabset",
		tabPanel("Correlation", 			
		isolate(selectInput("phen0", "Select Phenotype",choices=cn)),
		isolate(actionButton("go_alpha2", "Go!"))
	))
	appendTab(inputId = "tabset",
		tabPanel("Compare", 			
		isolate(selectInput("phen1", "Phenotype-1",choices=cn)),
		isolate(selectInput("phen2", "Phenotype-2",choices=cn)),
		isolate(selectInput("phen3", "Phenotype-3",choices=cn)),
		isolate(actionButton("go_alpha3", "Go!"))
	))
}

cmp_traits=function(v,my_trait1,my_trait2,my_trait3){
	print("INFO|cmp_traits")

	shinyjs::show("plot")
	shinyjs::hide("plot2")

	A=(colnames(v$transcriptomics))
	B=v$trait[,1]
	AB=B[B%in%A]
	ix=which(v$trait[,1]%in%AB)
	N=dim(v$transcriptomics)[1]
	SET_A=c(); SET_B=c(); SET_C=c()
	for(x in 1:N){
		if(my_trait1!="---"){
			iy=which(colnames(v$trait)==my_trait1)
			T1=as.numeric(v$trait[ix,iy])
			T2=as.numeric(unlist(v$transcriptomics[x,AB]))
			my_p=(cor.test(T1,T2))
			if(my_p$p.value<0.05){
				SET_A=c(SET_A,v$transcriptomics[x,1])
			}
		}
		if(my_trait2!="---"){
			iy=which(colnames(v$trait)==my_trait2)
			T1=as.numeric(v$trait[ix,iy])
			T2=as.numeric(unlist(v$transcriptomics[x,AB]))
			my_p=(cor.test(T1,T2))
			if(my_p$p.value<0.05){
				SET_B=c(SET_B,v$transcriptomics[x,1])
			}
		}
		if(my_trait3!="---"){
			iy=which(colnames(v$trait)==my_trait3)
			T1=as.numeric(v$trait[ix,iy])
			T2=as.numeric(unlist(v$transcriptomics[x,AB]))
			my_p=(cor.test(T1,T2))
			if(my_p$p.value<0.05){
				SET_C=c(SET_C,v$transcriptomics[x,1])
			}
		}

	}
	area1=SET_A[!(SET_A%in%SET_B)]
	area2=SET_B[!(SET_B%in%SET_A)]
	cross=SET_A[SET_A%in%SET_B]
	area1=length(area1)
	area2=length(area2)
	cross=length(cross)
	if(length(SET_C)>0){
		n12=SET_A[SET_A%in%SET_B]; n12=length(n12)
		n23=SET_B[SET_B%in%SET_C]; n23=length(n23)
		n13=SET_A[SET_A%in%SET_C]; n13=length(n13)
		n123=SET_A[SET_A%in%SET_B & SET_A%in%SET_C]; n123=length(n123)
		draw.triple.venn(area1=length(SET_A),area2=length(SET_B),area3=length(SET_C),n12,n23,n13,n123,
		 category        = c(my_trait1,my_trait2,my_trait3),
						fill            = c("blue","yellow","red"),
						lty             = "blank",
						cex             = 1,
						cat.cex         = 2,
						cat.dist        = -.5,
		)
	}else{
		draw.pairwise.venn(area1 = area1, area2=area2,cross.area=cross,
		 category        = c(my_trait1,my_trait2),
						fill            = c("blue","red"),
						lty             = "blank",
						cex             = 1,
						cat.cex         = 2,
						cat.dist        = -.5,
						ext.line.lwd    = 2,
					ext.line.lty    = "dashed")
	}
}

make_corr=function(v,my_trait,output){

	print("INFO|make_corr")

	shinyjs::show("plot")
	shinyjs::show("plot2")

	A=(colnames(v$transcriptomics))
	B=v$trait[,1]
	AB=B[B%in%A]
	ix=which(v$trait[,1]%in%AB)
	N=dim(v$transcriptomics)[1]
	X1=c(); X2=c()
	for(x in 1:N){
		iy=which(colnames(v$trait)==my_trait)
		T1=as.numeric(v$trait[ix,iy])
		T2=unlist(v$transcriptomics[x,AB])
		my_p=(cor.test(T1,T2))
		X1=c(X1,my_p$estimate)
		X2=c(X2,my_p$p.value)
		print(paste0("INFO|calc_cmp_transcriptomics_traits|",str(x),"|",str(N)))
	}
	L=list()
	L[[1]]=X1
	L[[2]]=p.adjust(X2,method="bonferroni")
	plot(L[[1]],L[[2]],log="y",xlab="correlation coefficient (r)",ylab="p.value")

	df=data.frame(L[[1]],L[[2]])
	df=df[order(df[,2]),]

	my_size=df[,2]

	n_all=length(my_size)
	n_5p_sign=length(my_size[my_size<0.05])
	n_1pm_sign=length(my_size[my_size<0.001])

	my_size=rep(1,length(df[,2]))
	my_size[df[,2]<0.05]=2
	my_size[df[,2]<0.001]=4
	
	L=list()
	L[[1]]=n_all
	L[[2]]=n_5p_sign
	L[[3]]=n_1pm_sign
			
	output$plot2=renderPlot({
		barplot(c(L[[1]],L[[2]],L[[3]]),col=c("green","#ff000044","#ff0000ff"),names=c("all genes","sign (p<0.05)","sign (p<0.001)"))
	})
	print(L)
	plot(df[,1],df[,2],pch=20,xlim=c(-1,1),cex=my_size,xlab="correlation coefficient (r)",ylim=c(0,1),ylab="p.value",main=my_trait)
	abline(h=0.05,lty=3)
	abline(h=0.001,lty=2)
	return(L)
}



