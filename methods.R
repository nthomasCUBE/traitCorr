library(d3heatmap)
library(gplots)
library(officer)
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
	cn2=c("---",v$transcriptomics[,1])
	appendTab(inputId = "tabset",
		tabPanel("Correlation with a trait", 			
		isolate(selectInput("phen0", "Select Phenotype",choices=cn)),
		radioButtons("corr_type", "Correlation coefficient:", c("Spearman" = "spearman","Pearson" = "pearson")),
		radioButtons("multiple_test_correction", "Multiple testing correction:", c("Benjamini-Hochberg (BH)"="BH","Bonferroni" = "bonferroni")),
		isolate(actionButton("go_alpha2", "Go!"))
	))
	appendTab(inputId = "tabset",
		tabPanel("Compare phenotypes", 			
		isolate(selectInput("phen1", "Phenotype-1",choices=cn)),
		isolate(selectInput("phen2", "Phenotype-2",choices=cn)),
		isolate(selectInput("phen3", "Phenotype-3",choices=cn)),
		radioButtons("corr_type2", "Correlation coefficient:", c("Spearman" = "spearman","Pearson" = "pearson")),
		radioButtons("multiple_test_correction2", "Multiple testing correction:", c("Benjamini-Hochberg (BH)"="BH","Bonferroni" = "bonferroni")),
		radioButtons("significance_level", "Significance level:", c("0.001"=0.001,"0.05" = 0.05)),
		isolate(actionButton("go_alpha3", "Go!"))
	))
	appendTab(inputId = "tabset",
		tabPanel("Regression analysis", 			
		isolate(selectInput("phen4", "Phenotype-1",choices=cn)),
		isolate(selectInput("gene1", "Gene-1",choices=cn2)),
		isolate(actionButton("go_alpha4", "Go!"))		
	))
	appendTab(inputId = "tabset",
		tabPanel("Download", 			
		downloadButton("download1","Sign corr (calculated so far, CSV)"),
		downloadButton("download2","Sign corr (calculated so far, DOCX)")
	))

}

cmp_traits=function(v,my_trait1,my_trait2,my_trait3){
	print("INFO|cmp_traits|start")

	shinyjs::show("plot")
	shinyjs::hide("plot2")

	mt_cor=v$multiple_test_correction2
	cor_type=v$corr_type2
	sign_level=as.double(v$significance_level)

	A=(colnames(v$transcriptomics))
	B=v$trait[,1]
	AB=B[B%in%A]
	ix=which(v$trait[,1]%in%AB)
	N=dim(v$transcriptomics)[1]
	SET_A=c(); SET_B=c(); SET_C=c()

	raw_p_val1=c()
	raw_p_val2=c()
	raw_p_val3=c()

	for(x in 1:N){
		print(paste0("INFO|cmp_traits|start",x,"|",N))
		print(sign_level)
		if(my_trait1!="---"){
			iy=which(colnames(v$trait)==my_trait1)
			T1=as.numeric(v$trait[ix,iy])
			T2=as.numeric(unlist(v$transcriptomics[x,AB]))
			my_p=(cor.test(T1,T2,method=cor_type))
			SET_A=c(SET_A,v$transcriptomics[x,1])
			raw_p_val1=c(raw_p_val1,my_p$p.value)
		}
		if(my_trait2!="---"){
			iy=which(colnames(v$trait)==my_trait2)
			T1=as.numeric(v$trait[ix,iy])
			T2=as.numeric(unlist(v$transcriptomics[x,AB]))
			my_p=(cor.test(T1,T2,method=cor_type))
			SET_B=c(SET_B,v$transcriptomics[x,1])
			raw_p_val2=c(raw_p_val2,my_p$p.value)
		}
		if(my_trait3!="---"){
			iy=which(colnames(v$trait)==my_trait3)
			T1=as.numeric(v$trait[ix,iy])
			T2=as.numeric(unlist(v$transcriptomics[x,AB]))
			my_p=(cor.test(T1,T2,method=cor_type))
			SET_C=c(SET_C,v$transcriptomics[x,1])
			raw_p_val3=c(raw_p_val3,my_p$p.value)
		}

	}
	print(summary(raw_p_val1))
	print(summary(raw_p_val2))
	print(summary(raw_p_val3))
	
	raw_p_val1=p.adjust(raw_p_val1,mt_cor)
	SET_A=SET_A[raw_p_val1<=sign_level]
	raw_p_val2=p.adjust(raw_p_val2,mt_cor)
	SET_B=SET_B[raw_p_val2<=sign_level]
	raw_p_val3=p.adjust(raw_p_val3,mt_cor)
	SET_C=SET_C[raw_p_val3<=sign_level]
	
	print(SET_A)
	print(SET_B)
	print(SET_C)
	
	print(paste0("SET_A","|",length(SET_A)))
	print(paste0("SET_B","|",length(SET_B)))
	print(paste0("SET_C","|",length(SET_C)))
	
	area1=SET_A[!(SET_A%in%SET_B)]
	area2=SET_B[!(SET_B%in%SET_A)]
	cross=SET_A[SET_A%in%SET_B]
	area1=length(area1)
	area2=length(area2)
	cross=length(cross)
	
	if((length(SET_A)+length(SET_B)+length(SET_C))==0){
		shinyalert("INFO", "No genes remained using the current filters!", type = "info")
	}else if(length(SET_C)>0){
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

regr_analysis=function(v, my_trait1, my_gene1){
	print("INFO|regr_analysis")
	ix=which(my_gene1==v$transcriptomics[,1])

	data1=v$transcriptomics
	data2=v$trait
	
	A=colnames(data1)
	B=data2[,1]
	AB=A[A%in%B]
	print(length(AB))

	ix=which(data1[,1]==my_gene1)
	U1=(as.numeric(data1[ix,AB]))

	rownames(data2)=data2[,1]
	U2=(data2[AB,my_trait1])

	ct=cor.test(U1,U2,method="spearman")
	
	plot(U1,U2,xlab=my_gene1,ylab=my_trait1)
	abline(lm(U2~U1))
	legend("bottom",paste0("pvalue=",ct$p.value))
	legend("top",paste0("r=",ct$estimate))

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
	X1=c(); X2=c(); X3=c(); X4=c();
	for(x in 1:N){
		iy=which(colnames(v$trait)==my_trait)
		T1=as.numeric(v$trait[ix,iy])
		T2=unlist(v$transcriptomics[x,AB])
		my_p=(cor.test(T1,T2,method=v$corr_type))
		X1=c(X1,my_p$estimate)
		X2=c(X2,my_p$p.value)
		X3=c(X3,v$transcriptomics[x,1])
		X4=c(X4,my_trait)
	}
	L=list()
	X2=p.adjust(X2,method=v$multiple_test_correction)
	L[[1]]=X1
	L[[2]]=X2
	L[[3]]=X3
	L[[4]]=X4
	plot(L[[1]],L[[2]],log="y",xlab="correlation coefficient (r)",ylab="p.value")
	new_entry=data.frame()
	for(x in 1:length(X1)){
		if(X2[x]<0.05){
			new_entry=rbind(new_entry,c(X1[x],X2[x],X3[x],X4[x]))
		}
	}

	df_all=data.frame(L[[1]],L[[2]],L[[3]],L[[4]])
	df_all=subset(df_all,df_all[,2]<0.001)
	gene_selected=(df_all[,3])
	u_m=unique(v$module[,2])
	m_arr=c()
	m_arr_1=c()
	m_arr_2=c()
	for(x in 1:length(u_m)){
		g_a=subset(v$module,v$module[,2]==u_m[x])[,1]
		g_s=g_a[g_a%in%gene_selected]
		m_r=100*length(g_s)/length(g_a)
		m_arr=c(m_arr,m_r)
		m_arr_1=c(m_arr_1,length(g_s))
		m_arr_2=c(m_arr_2,length(g_a))
	}
	
	
	output$plot3=renderPlot({
		#barplot(m_arr,names=u_m)
		
		A=m_arr_1
		B=m_arr_2
		
		m3=c()
		for(x in 1:length(A)){
			x1=A[x]
			x2=sum(A)-x1
			x3=B[x]
			x4=sum(B)-x3
			print(c(x1,x2,x3,x4))
			m3=c(m3,(phyper(A[x],B[x],sum(B),sum(A),lower.tail=FALSE)))
		}
		m3=-log2(m3)
		m3[m3>10]=10
		print(m3)
		
		print(m3)
		plot(m_arr_1,m_arr_2,cex=m3,pch=20)
		text(m_arr_1,m_arr_2,u_m,cex=2.5)
		
		print(m_arr_1)
		print(m_arr_2)
	})
	

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
	plot(df[,1],df[,2],pch=20,xlim=c(-1,1),cex=my_size,xlab="correlation coefficient (r)",ylim=c(0,1),ylab="p.value",main=my_trait)
	abline(h=0.05,lty=3)
	abline(h=0.001,lty=2)

	v$df_output=new_entry

	return(L)
}



