
library(readr) # for files importing
library(miRBaseConverter) # for miRNA'names conversion
library(dplyr) # for data manipulation

############### processing MTI data ############### 
# Set as working directory the folder containing your prediction/interaction data. Here
# we'll use experimentally validated miRNA-target interaction (MTI) 
# data from TarBase downloaded from https://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Fdownloaddataform
# you have to fill the form to download it (they will send you a link to your academic institution email)

setwd("D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/predictions")

# import interaction data
tarbase <- as.data.frame(read_delim("tarbase.txt", "\t", escape_double = TRUE, 
                                    trim_ws = TRUE))

# Restrict interactions filtering to Homo sapiens and DIRECT validation type 
tarbase <- tarbase[tarbase$species=="Homo sapiens" & tarbase$direct_indirect=="DIRECT",]

# Reorder and keep only columns 2 and 3
tarbase = tarbase[,c(3,2)]
# rename it so we can merge it later
names(tarbase)[1] <- "a"
names(tarbase)[2] <- "b"

############### processing MTI data ############### 
# Set as working directory the folder containing your prediction/interaction data. Here
# we'll use experimentally validated miRNA-target interaction (MTI) 
# data from miRTarBase downloaded from http://mirtarbase.cuhk.edu.cn/php/download.php
setwd("D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/predictions")
# import interaction data
mirtarbase <- as.data.frame(read_delim("mirtarbase.txt", "\t", escape_double = TRUE, 
                    trim_ws = TRUE))

# Reorder and keep only columns 2 and 3
mirtarbase = mirtarbase[,c(2,4)]
# rename it so we can merge it later
names(mirtarbase)[1] <- "a"
names(mirtarbase)[2] <- "b"

# merge data frames by rows
mtb = rbind(mirtarbase, tarbase)

# adjust data frames names again
names(mtb)[names(mtb) == "a"] <- "miRNA"
names(mtb)[names(mtb) == "b"] <- "Genes"

# create a vector with miRNA's names
miRNANames <- mtb$miRNA

# convert the its mature miRNA name to precursor so we can cross these interactions
# with expression data and infer some regulation
precursors <- miRNA_MatureToPrecursor(miRNANames)

# merge the precursors dataframe with the old data frame first imported 
# excluding the first column (which is also present in precursor data.frame)
# also excluind 5-8 because it's not useful for us now
converted_miRNAs <- cbind (precursors, mtb[2])

# filter to keep only miRNAs that physically interact with PAWR (or your gene of interest)
PAWR <- unique(converted_miRNAs[converted_miRNAs$Genes=="PAWR",])
### consider using the code below if you studying many genes (specifying genes names in 
# a vector class object or cocatenated with c() as shown below):

# miRNA_filt <- 
  # mirNA_converted[mirNA_converted$`Target Gene` %in% c(), ]



# create a vector only with miRNAs precursor names
precursor_mirs <- as.vector(unique(PAWR[,2]))

# Now make another vector of miRNAs by its mature name 
mature_mirs <- as.vector(unique(PAWR[,1]))

# remove unecessary objects
rm(precursors, miRNANames, PAWR)



############### crossing miRNA expression data and its interactions  ############### 
# Now we are going to cross the interaction information filtered before with expression data

# Set the folder containing the expression file that we downloaded from TCGA and analyzed 
# for DEGs using edgeR
setwd('D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/resultados/mirna_estratificacao')

# import expression data file
edger <- read_delim("toptags.txt", "\t", escape_double = TRUE, 
                         trim_ws = TRUE)

# Filter miRNAs that interact with PAWR
pawr_reg <- edger[edger$X1 %in% precursor_mirs, ]



# Create a data frame of miRNAs that interact with PAWR and are up/down-regulated
# with it's respective expression values
up <- pawr_reg[pawr_reg$logFC > 1.5 & pawr_reg$FDR < 0.05 ,]
down <- pawr_reg[pawr_reg$logFC < -1.5 & pawr_reg$FDR < 0.05 ,]

# Filter miRNA-mRNA interactions (without expression values) to build a interaction network
# of miRNA that target PAWR and other genes
# Notice that here not only miRNA-PAWR interactions are included
filt_up <- unique(converted_miRNAs[converted_miRNAs$Precursor %in% up$X1, ])
filt_down <- unique(converted_miRNAs[converted_miRNAs$Precursor %in% down$X1, ])
# Actually, we are going to use that to build a network on Enrichment Map cytoscape tool


# Keep miRNAs that really interact with PAWR (according to miRTarBase data) 
# by its mature names. Notice that we have to do that because of mature > precursor step
up_mat <- filt_up[filt_up$OriginalName %in% mature_mirs, ]
down_mat <- filt_down[filt_down$OriginalName %in% mature_mirs, ]

rm (pawr_reg, precursor_mirs, converted_miRNAs, filt_down, filt_up, mature_mirs)





############### Restrict data to some genes of interest ############### 
# if you want to build a interaction network of specific genes, do this:
setwd("D:/Faculdade/pibit/PIBIT/Sergio/Continuacao_EM_PAWR/ENRICHMENT-MAP_MIRNA-FUNCOES-EXP/rede_mirna-genes_sem-funcao-reg-gene-exp/")
interest <- as.list(read_delim("genes_de_interesse.txt", "\t", escape_double = TRUE, 
                         trim_ws = TRUE))

View(interest) #notice that the list has a header named 'Genes'

up_int <- UP_mat[UP_mat$`Target Gene` %in% interest$Genes, ]

up_int <- UP_mat[UP_mat$`Target Gene` %in% interest$Genes, ]




############### Save results ############### 
setwd("D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/resultados/mirna_estratificacao")

dflist <- list(down_mat, up_mat)
names(dflist) <- c("down_mat", "up_mat")
dflist_names <- names(dflist)

for(i in 1:length(dflist)){
  write.table(dflist[[i]],
              paste0(dflist_names[i], ".txt"), 
              sep = "\t",
              row.names = TRUE, 
              col.names = NA, 
              qmethod = c("escape", "double"))
}

rm (dflist, dflist_names, i)

