#######PCA By Tumor Sample
Data = read.table("GSE115719_Tumor_raw_count.txt", header=TRUE) ###Input Data
GeneID = Data[,2]## Save GeneIDs
Data = Data[3:22]  ###Grab only the tumor samples
##Sum up technical replicates
colnames(Data)<- c()## Remove header
M=data.frame(Data[,1]+Data[,2])###Initilize with 1st calc
for(i in seq(from=3, to=19, by=2)){
	x = data.frame(Data[,i]+Data[,i+1])
	M <- data.frame(M,x)
}
Data = M
library("DESeq2")
dds <- DESeqDataSetFromMatrix(Data = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)

Data <- counts(dds, normalized=TRUE)
colnames(Data)<- c()## Remove header
Data = as.matrix(Data) ###make into matrix
data.pca<-prcomp(Data,center=TRUE,scale=TRUE) ###perform PCA
plot(data.pca, type="l")### plot varaince of PC
### Variance 
pr.var=data.pca$sdev ^2 
pr.var 
pve=pr.var/sum(pr.var)
pve


sample=c(1,2,3,4,5,6,7,8,9,10) ###Set sample numbers
s = data.frame(sample,data.pca$rotation) ###append sample numbers with PCA loadings. 

###Plotting

##Color Pallet
colors <- c("#999999", "#E69F00", "#56B4E9", "#0000CD", "#3CB371", "#C0FF3E", "#FFFF00", "#EE4000", "#1E1E1E", "#A2CD5A") 
colors <- colors[as.numeric(s$sample)]

##Plots
#with(s,plot(PC1,PC2,col=colors,pch=19))  
#with(s,plot(PC1,PC3,col=colors,pch=19)) 
#with(s,plot(PC2,PC3,col=colors,pch=19)) 
with(s,plot(PC4,PC3,col=colors,pch=19)) 

##Legend
legend("topleft", inset=.05,      # location and inset
    bty="n", cex=.5,              # suppress legend box, shrink text 50%
    title="Tumor Samples",
    c(" Tumor 1", " Tumor 2", " Tumor 3", " Tumor 4 ", " Tumor 5 ", " Tumor 6 ", " Tumor 7 ", " Tumor 8 ", " Tumor 9 ", " Tumor 10 "), fill=c("#999999", "#E69F00", 
"#56B4E9", "#0000CD", "#3CB371", "#C0FF3E", 
"#FFFF00", "#EE4000", "#1E1E1E", "#A2CD5A"))

biplot(data.pca , scale=0)##biplot 



#######PCA By Gene
s = as.matrix(t(Data)) ##Tranpose Data
GeneID = (t(GeneID)) ##Tranpose GeneID
data.pca = prcomp(s)
plot(data.pca, type="l")### plot varaince of PC
### Variance 
pr.var=data.pca$sdev ^2 
pr.var 
pve=pr.var/sum(pr.var)
pve
g = data.frame(GeneID,data.pca$rotation) ###append GeneID with PCA loadings. 
##Plots
with(g,plot(PC1,PC2))  
with(g,plot(PC1,PC3)) 
with(g,plot(PC2,PC3)) 
with(g,plot(PC4,PC3)) 




