library(d3heatmap)
library(scales)
library(xlsx)
library(gplots)

options(stringsAsFactors=FALSE)

add_transcriptomics=function(){
	print(c("INFO|add_transcriptomcis"))
}

add_trait_information=function(){
	print(c("INFO|add_trait_information"))
}

calc_cmp_transcriptomics_traits=function(v){
	print(c("INFO|transcriptomics VS traits"))
	
	cn=colnames(v$trait)
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

cmp_traits=function(v,my_trait1,my_trait2){
	A=(colnames(v$transcriptomics))
	B=v$trait[,1]
	AB=B[B%in%A]
	ix=which(v$trait[,1]%in%AB)
	N=dim(v$transcriptomics)[1]
	SET_A=c(); SET_B=c()
	for(x in 1:N){
		iy=which(colnames(v$trait)==my_trait1)
		T1=as.numeric(v$trait[ix,iy])
		T2=unlist(v$transcriptomics[x,AB])
		my_p=(cor.test(T1,T2))
		if(my_p<0.05){
			SET_A=c(SET_A,v$transcriptomics[x,1])
		}
		iy=which(colnames(v$trait)==my_trait2)
		T1=as.numeric(v$trait[ix,iy])
		T2=unlist(v$transcriptomics[x,AB])
		my_p=(cor.test(T1,T2))
		if(my_p<0.05){
			SET_B=c(SET_B,v$transcriptomics[x,1])
		}
	}
	print(SET_A)
	print(SET_B)
}

make_corr=function(v,my_trait){
	A=(colnames(v$transcriptomics))
	B=v$trait[,1]
	AB=B[B%in%A]
	ix=which(v$trait[,1]%in%AB)
	N=dim(v$transcriptomics)[1]
	X1=c()
	X2=c()
	for(x in 1:N){
		iy=which(colnames(v$trait)==my_trait)
		print(paste0("my_trait::",my_trait))
		print(paste0("tratis::",colnames(v$trait)))
		T1=as.numeric(v$trait[ix,iy])
		T2=unlist(v$transcriptomics[x,AB])
		print(T1)
		print(T2)
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
	print(df[,2])
	my_size=rep(1,length(df[,2]))
	my_size[df[,2]<0.05]=2
	my_size[df[,2]<0.001]=4
	print(my_size)
	plot(df[,1],df[,2],pch=20,cex=my_size,xlab="correlation coefficient (r)",ylim=c(0,1),ylab="p.value",main=my_trait)
	abline(h=0.05,lty=3)
	abline(h=0.001,lty=2)
}



