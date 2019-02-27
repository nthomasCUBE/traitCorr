d1=read.csv("data/FemaleLiver-Data/ClinicalTraits.csv",sep=",",header=T,row.names=1)
write.table(d1,"Traits_modified.txt",sep="\t",row.names=FALSE,quote=FALSE)

d1=read.csv("data/FemaleLiver-Data/LiverFemale3600.csv",sep=",",header=T)
d1=cbind("Gene"=rn,d1[,9:dim(d1)[2]])
write.table(d1,"Expression_modified.txt",sep="\t",row.names=FALSE,quote=FALSE)

