#setwd("path/to/your/dataframe/file) #MacOS
#setwd("path\\to\\your\\dataframe\\file) #Windows

### IF USED ON WINDOWS, CHANGE ALL / BY \\ FOR PATHFILES
setwd("/Users/adam-nicolaspelletier/Desktop/rQTL_script") # SEt up your working directry, i.e. the location where your csv fiel is located for your analysis
filename = "SNPADAM.csv"

redoperms = "NO" ## Change to YES if you wish to recompute the permutations

library(qtl)

if (dir.exists("Figures") == FALSE) {
  dir.create("Figures")
}
if (dir.exists("Output_dataframes") == FALSE) {
  dir.create("Output_dataframes")
}

read.cross(format="csvr", file=filename, na.strings=c("-","NC"), genotypes=c("0", "1", "2", "D", "C"))->nb10
colnames(nb10$pheno) -> phenotypenames

calc.genoprob(nb10, step=4) -> nb10


sc1.test <- scanone(nb10, pheno.col=phenotypenames[1], model="normal", method="mr")
perms <- scanone(nb10, method=”mr”, model, method="normal", normal, pheno.col=1, n.perm=1000)
perm.mr[[i]] = summary(perms, alpha=c(0.65, 0.05))

sc1.mr <- list()
sc1.hk <- list()
sc2.mr <- list()
sc2.hk <- list()



##Scanone 
if (redoperms == "YES") {
  perms <- scanone(nb10, method="mr", model= "normal", pheno.col=2, n.perm=10000)
  perm.sig = summary(perms, alpha=c(0.05))[[1]]
  perm.sug = summary(perms, alpha=c(0.65))[[1]]
  perms2 <- scantwo(nb10, method="mr", model= "normal", pheno.col=2, n.perm=1000)
  perm2.sig = summary(perms2, alpha=c(0.05))[[1]] ## 1 for full sc2 lodscore thresholds. view cantwo documentation for other lodscore options
  perm2.sug = summary(perms2, alpha=c(0.65))[[1]]
}


for(i in 1:length(phenotypenames)){
  sc1.mr[[i]] = scanone(nb10, pheno.col=i, model="normal", method="mr")
  sc1.hk[[i]] = scanone(nb10, pheno.col=i, model="normal", method="hk")
}


## Scantwo
for(i in 1:length(phenotypenames)){
  sc2.mr[[i]] = scantwo(nb10, pheno.col=i, model="normal", method="mr")
  sc2.hk[[i]] = scantwo(nb10, pheno.col=i, model="normal", method="hk")
}

figformat = "svg" #svg, or tiff. can also add other formats such as png or pdf in the code if needed

################ Output for sc1.mr#########################################
sc1.mr.header = paste("Significant LOD Score:", perm.sig,"\nSuggestive LOD Score:", perm.sug,"\n", sep="")
cat(sc1.mr.header,file= paste("Output_dataframes/sc1.mr.txt"))

for(i in 1:(length(phenotypenames)-1)){
  filename.mr = paste(i,"_",phenotypenames[i],".sc1.mr", sep="")
  if (figformat=="svg"){
    svg(file=paste("Figures/",filename.mr,".svg", sep=""), width = 8, height = 7)
  } else if (figformat =="tiff"){
      tiff(file=paste("Figures/",filename.mr,".tiff", sep=""), width = 480, height = 480,
             units = "px")
  }
  
  plot(sc1.mr[[i]], chr=-20)
  abline(h=perm.sig,lty =2)
  abline(h=perm.sug,lty =3)
  dev.off()
  
  if (i < 2) {
    sc1.mr.total = sc1.mr[[i]]
    rownames(sc1.mr.total) -> sc1.mr.total$names
    names(sc1.mr.total)[names(sc1.mr.total)=="lod"] <- phenotypenames[i] 
  } else if (i >= 2) {
    df = sc1.mr[[i]]
    rownames(df) -> df$names
    names(df)[names(df)=="lod"] <- phenotypenames[i]
    df$pos <- NULL
    df$chr <- NULL
    merge(sc1.mr.total,df,"names") -> sc1.mr.total
  }
}

sc1.mr.totalf <- sc1.mr.total[order(sc1.mr.total$chr,sc1.mr.total$pos),]
write.table(sc1.mr.totalf, file=paste("Output_dataframes/sc1.mr.txt"), append=TRUE,
              sep="\t",col.names=NA)

################ Output for sc1.hk#########################################
sc1.hk.header = paste("Significant LOD Score:", perm.sig,"\nSuggestive LOD Score:", perm.sug,"\n", sep="")
cat(sc1.hk.header,file= paste("Output_dataframes/sc1.hk.txt"))

for(i in 1:(length(phenotypenames)-1)){
  filename.hk = paste(i,"_",phenotypenames[i],".sc1.hk", sep="")
  if (figformat=="svg"){
    svg(file=paste("Figures/",filename.hk,".svg", sep=""), width = 8, height = 7)
  } else if (figformat =="tiff"){
    tiff(file=paste("Figures/",filename.hk,".tiff", sep=""), width = 480, height = 480,
         units = "px")
  }
  
  plot(sc2.hk[[i]], chr=-20)
  abline(h=perm.sig,lty =2)
  abline(h=perm.sug,lty =3)
  dev.off()
  
  if (i < 2) {
    sc1.hk.total = sc1.hk[[i]]
    rownames(sc1.hk.total) -> sc1.hk.total$names
    names(sc1.hk.total)[names(sc1.hk.total)=="lod"] <- phenotypenames[i] 
  } else if (i >= 2) {
    df = sc1.hk[[i]]
    rownames(df) -> df$names
    names(df)[names(df)=="lod"] <- phenotypenames[i]
    df$pos <- NULL
    df$chr <- NULL
    merge(sc1.hk.total,df,"names") -> sc1.hk.total
  }
}

sc1.hk.totalf <- sc1.hk.total[order(sc1.hk.total$chr,sc1.hk.total$pos),]
write.table(sc1.hk.totalf, file=paste("Output_dataframes/sc1.hk.txt"), append=TRUE,
            sep="\t",col.names=NA)

################ Output for sc2.mr#########################################
sc2.mr.header = paste("Significant LOD Score:", perm2.sig,"\nSuggestive LOD Score:", perm2.sug,"\n", sep="")
cat(sc2.mr.header,file= paste("Output_dataframes/sc2.mr.txt"))

for(i in 1:(length(phenotypenames)-1)){
  filename.2.mr = paste(i,"_",phenotypenames[i],".sc2.mr", sep="")
  if (figformat=="svg"){
    svg(file=paste("Figures/",filename.2.mr,".svg", sep=""), width = 8, height = 7)
  } else if (figformat =="tiff"){
    tiff(file=paste("Figures/",filename.2.mr,".tiff", sep=""), width = 480, height = 480,
         units = "px")
  }
  
  plot(sc2.mr[[i]], chr=-20)
  dev.off()
  
  if (i < 2) {
    sc2.mr.total = summary(sc2.mr[[i]])
    rownames(sc2.mr.total) -> sc2.mr.total$names
    names(sc2.mr.total)[names(sc2.mr.total)=="lod.full"] <- phenotypenames[i] 

  } else if (i >= 2) {
    df = summary(sc2.mr[[i]])
    rownames(df) -> df$names
    names(df)[names(df)=="lod.full"] <- phenotypenames[i]
    df$chr1<- NULL
    df$chr2 <- NULL
    df$pos1f <- NULL
    df$pos2f <- NULL
    df$lod.fv1 <- NULL
    df$lod.int <- NULL
    df$pos1a <- NULL
    df$pos2a <- NULL
    df$lod.add <- NULL
    df$lod.av1 <- NULL
    merge(sc2.mr.total,df,"names") -> sc2.mr.total
  }
}

sc2.mr.total$pos1f <- NULL
sc2.mr.total$pos2f <- NULL
sc2.mr.total$lod.fv1 <- NULL
sc2.mr.total$lod.int <- NULL
sc2.mr.total$pos1a <- NULL
sc2.mr.total$pos2a <- NULL
sc2.mr.total$lod.add <- NULL
sc2.mr.total$lod.av1 <- NULL

sc2.mr.totalf <- sc2.mr.total[order(sc2.mr.total$chr1,sc2.mr.total$chr2),]
write.table(sc2.mr.totalf, file=paste("Output_dataframes/sc2.mr.txt"), append=TRUE,
            sep="\t",col.names=NA)


################ Output for sc2.hk#########################################
sc2.hk.header = paste("Significant LOD Score:", perm2.sig,"\nSuggestive LOD Score:", perm2.sug,"\n", sep="")
cat(sc2.hk.header,file= paste("Output_dataframes/sc2.hk.txt"))


cat(,file="outfile.txt",sep="\n")
for(i in 1:(length(phenotypenames)-1)){
  filename.2.hk = paste(i,"_",phenotypenames[i],".sc2.hk", sep="")
  if (figformat=="svg"){
    svg(file=paste("Figures/",filename.2.hk,".svg", sep=""), width = 8, height = 7)
  } else if (figformat =="tiff"){
    tiff(file=paste("Figures/",filename.2.hk,".tiff", sep=""), width = 480, height = 480,
         units = "px")
  }
  
  plot(sc2.hk[[i]], chr=-20, upper=NULL)
  dev.off()
  
  if (i < 2) {
    sc2.hk.total = summary(sc2.hk[[i]])
    rownames(sc2.hk.total) -> sc2.hk.total$names
    names(sc2.hk.total)[names(sc2.hk.total)=="lod.full"] <- phenotypenames[i] 
    
  } else if (i >= 2) {
    df = summary(sc2.hk[[i]])
    rownames(df) -> df$names
    names(df)[names(df)=="lod.full"] <- phenotypenames[i]
    df$chr1<- NULL
    df$chr2 <- NULL
    df$pos1f <- NULL
    df$pos2f <- NULL
    df$lod.fv1 <- NULL
    df$lod.int <- NULL
    df$pos1a <- NULL
    df$pos2a <- NULL
    df$lod.add <- NULL
    df$lod.av1 <- NULL
    merge(sc2.hk.total,df,"names") -> sc2.hk.total
  }
}

sc2.hk.total$pos1f <- NULL
sc2.hk.total$pos2f <- NULL
sc2.hk.total$lod.fv1 <- NULL
sc2.hk.total$lod.int <- NULL
sc2.hk.total$pos1a <- NULL
sc2.hk.total$pos2a <- NULL
sc2.hk.total$lod.add <- NULL
sc2.hk.total$lod.av1 <- NULL

sc2.hk.totalf <- sc2.hk.total[order(sc2.hk.total$chr1,sc2.hk.total$chr2),]
write.table(sc2.hk.totalf, file=paste("Output_dataframes/sc2.hk.txt"), append=TRUE,
            sep="\t",col.names=NA)









