d1=read.csv("data/FemaleLiver-Data/ClinicalTraits.csv",sep=",",header=T,row.names=1)
write.table(d1,"Traits_modified.txt",sep="\t",row.names=FALSE,quote=FALSE)

