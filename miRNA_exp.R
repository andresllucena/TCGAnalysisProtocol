setwd("D:/Faculdade/pibit/PIBIT/Sergio/estratificacao/sample10")

############### Do edgeR DEG analysis ############### 

# Install and import edgeR
install.packages('edgeR')
library(edgeR) 

x = list()
et = list()
tt = list()

for (i in 1:10){
  setwd(file.path('D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/resultados/mirna_estratificacao/', SampleNames[[i]]))
  # Merge the two groups of tissue samples (normal/cancer) into one data.frame
  y <- cbind(Normal_miRNA, StrSample[[i]])  

  # Now define the group of each sample
  group <- c()
  group_num <- c(1,2)
  group[1:46]<- unique(group_num)[1] 
  group[47:196]<- unique(group_num)[2] 
  
  
  # Create a DGElist object
  bb <- DGEList(counts=y[1:196], group=group)
  #x[["counts"]] 
  x[[i]] <- bb

  
  # Filter genes low expressed in all samples
  z <- filterByExpr(x[[i]], min.count = 0.1, min.total.count = 0.1)
  x[[i]] <- x[[i]][z, , keep.lib.sizes=FALSE]
  nrow(x[[i]]) # see how many genes it still has
  table(z) # see how many genes were FALSE and TRUE for the lib size default conditions
  x[[i]]$samples$lib.size <- colSums(x[[i]]$counts) # recompute the library sizes

  ############### Normalization and exploratory analysis ############### 
  
  x[[i]] <- calcNormFactors(x[[i]])
  
  # Create a Multidimensional Scaling Plot (plotMDS)
  colors <- c("green", "red") # choose colors according to the sample group order used in y. Here, Adenocarcinoma samples are red.
  points <- rep(c(20,18),2) # choose points style according to the same order
  png(filename= paste0(i,'plotMDS.png'))
  plotMDS <- plotMDS(x[[i]], top = 500, pch=points[group], col=colors[group], gene.selection = "pairwise", main="Similaridade de expressão entre amostras normais e tumorais (miRNA)") #da os nomes de acordo com os grupos
  group_name <- c("Tecido Normal", "Adenocarcinoma") # specify names to be used in legend for each sample group
  legend("topright", legend=group_name, 
         pch=points, col=colors, ncol=1, cex = 0.8)
  dev.off()
  
  rm(colors, group_name, plotMDS, points, z, y)

  
  # Biological Coefficient of Variation plot
  x[[i]] <- estimateCommonDisp(x[[i]]) 
  x[[i]] <- estimateTagwiseDisp(x[[i]])
  x[[i]] <- estimateTrendedDisp(x[[i]])
  # x$common.dispersion    ### Delete the # if you want to print it.
  # View(x[["tagwise.dispersion"]])
  # View(x[["trended.dispersion"]])
  png(filename= paste0(i,'plotBCV.png'))
  plotBCV(x[[i]])
  dev.off()
  
  png(filename= paste0(i,'plotMeanVar.png'))
  plotMeanVar(x[[i]], show.tagwise.vars=TRUE, NBline=TRUE)
  dev.off()
  # these plots give insights about quality control of the data to be analyzed
 
  
  ############### DE analysis  and it's related plots ############### 
  
  et[[i]] <- exactTest(x[[i]])
  topTags(et[[i]]) #show the top DE genes 
  
  # Create a TopTags list-type object with all DE analysis values and info about the method used. There is a table with all values.
  tt[[i]] = topTags(et[[i]], n=nrow(x[[i]]))
  head(tt[[i]]$table)
  
  # want to save a table from this file so you can explore results outside R?
  write.table(tt[[i]]$table, 
              file = paste0(i, "toptags.txt", sep = ""), sep = "\t",
              row.names = TRUE, col.names = NA, 
              qmethod = c("escape", "double"))
  

  ### Create a MD (mean-difference) plot of your expression data
  # See how many genes are down, up or do not have significative p-value
  summary(decideTests(et[[i]], adjust.method = "fdr", lfc = 2))
  
  # set adjustment method and fold change cut-off 
  is.de <- decideTestsDGE(et[[i]], adjust.method = "fdr", lfc = 2)
  
  # Choose dots appearance and some of its legend settings
  group_name <- c("Up-regulated (> 2)", "Down-regulated (< -2)")
  colors <- c("red", "blue")
  points <- 1
  
  png(filename = paste0(i,'plotMD.png'))
  # Make the graph
  plotMD(et[[i]], hl.pch=1, hl.col = colors, status = is.de, main = "Normal x Adenocarcinoma (miRNA)")
  # Overlap the cut-off lines (The green lines indicate 1-fold changes.)
  abline(h=c(-2, 2), col="green") 
  # Overlap a customized legend if you want to
  legend("topright", legend=group_name, 
         pch=1, col=colors, ncol=1, cex = 0.8)
  dev.off()
  
  setwd('D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/resultados/mirna_estratificacao')
}

rm(colors, group, group_name, group_num, points, bb)

############### Take means of toptags of all samples ############### 
# first, merge all data-frames of 'tt' together
zz = tt[[1]]$table
colnames(zz) <- paste(colnames(zz), 'a1', sep = "_") # change its names so 'merge' won't give warnings about duplicated columns
zz = tibble::rownames_to_column(zz, "Names") # convert row names to a column to avoid errors in merge function

# now, do this to all the data frames in tt
for(i in 2:10){
  a = tt[[i]]$table
  b = paste0('a', i, sep = '')
  colnames(a) <- paste(colnames(a), b, sep = "_")
  a = tibble::rownames_to_column(a, "Names") #
  zz = merge(zz, a, by = 1, all = TRUE)
}
rm(a, b, i, is.de, x)

# turn #1 column back to its place as rownames
zz = tibble::column_to_rownames(zz, "Names")

# make columns seem duplicated so we can do rowmeans on it later
names(zz) = gsub(pattern = c(as.character("_a.*")), replacement = "", x = names(zz)) #tirar tudo que tiver depois de 11A

# Now, take mean values of duplicated columns
zz <- as.data.frame( # sapply returns a list here, so we convert it to a data.frame
  sapply(unique(names(zz)), 
         function(col) rowMeans(zz[names(zz) == col]) 
  )
)

# save the top tags table in in your computer for further analysis of prediction and GO
write.table(zz, file = "toptags.txt", sep = "\t",
            row.names = TRUE, col.names = NA, 
            qmethod = c("escape", "double"))
