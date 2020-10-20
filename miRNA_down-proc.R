############### Import packages for initial data download and organization ###############

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("TCGAWorkflow")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks", force = TRUE)

browseVignettes("TCGAWorkflow")

# Download packages
pkgs = c('TCGAWorkflowData', 'DT', 'TCGAbiolinks', 
         'SummarizedExperiment', 'dplyr', 'readr', 'RMySQL')
install.packages(pkgs) # before go to shell and type '$ sudo apt-get install libmariadb-client-lgpl-dev'

# Load downloaded packages 
lapply(pkgs, require, character.only = TRUE)

# Delete vector with the name of the packages
rm(pkgs)


############### Downloading TCGA Lung cancer miRNA expression files ############### 
setwd("D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/resultados/mirna_estratificacao")

TCGAbiolinks:::getGDCprojects()$project_id  # gives a list of GDC projects IDs according to cancer/sample/data source types
TCGAbiolinks:::getProjectSummary("TCGA-LUAD") # gives a list with possible data categories of files of lung adenocarcinoma (LUAD) patients from TCGA that can be downloaded

queryLUAD.exp <- GDCquery(project = "TCGA-LUAD", # choose a TCGA project of interest
                          data.category = "Transcriptome Profiling", # (Biospecimen, Clinical, CNV, DNA methylation, Sequencing Reads, SNV, transcriptome profiling, Gene expression, Protein expression, Raw microarray data, Raw sequencing data)
                          data.type = "miRNA Expression Quantification", # data type to filter the files to download (Gene Expression Quantification, Isoform Expression Quantification, miRNA Expression Quantification, Splice Junction Quantification)
                          workflow.type = "BCGSC miRNA Profiling", # specify what type of counts you want (HTSeq-Counts, HTSeq-FPKM, HTSeq-FPKM-UQ, STAR-Counts)
                          experimental.strategy = "miRNA-Seq", #Filter to experimental strategy (WXS, RNA-Seq, miRNA-Seq, Genotyping Array, DNA-Seq, Methylation array, Protein expression array, CGH array, VALIDATION, Gene expression array, WGS, MSI-Mono-Dinucleotide Assay, miRNA expression array, Mixed strategies, AMPLICON, Exon array, Total RNA-Seq, Capillary sequencing, Bisulfite-Seq)
                          sample.type = c("Primary Tumor")) # you can find sample type complete list on https://www.bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/query.html#sampletype_options

GDCdownload(queryLUAD.exp) # then download the data to your computer (can be found in the active working directory)
luad.exp <- GDCprepare(query = queryLUAD.exp, save = TRUE, save.filename = "luadExp.rda") # and import it to your R enviroment
# The data downloaded here is from TCGA database. The LUAD project has data of lung adenocarcinoma tissue samples. We downloaded gene expression RNA-seq data from 519 samples (Primary solid Tumor) of lung adenocarcinoma.

#### Adjust columns names exchanging hyphens to underlines
# Assign with 'xxx' the columns that will be removed (we'll only keep with just raw counts columns)
names(luad.exp) = gsub(x = names(luad.exp), pattern = "\\-", replacement = "_") # 
names(luad.exp) = gsub(pattern = c(as.character("*.eads_per_million|*.ross_mapped")), 
                                replacement = "xxx", x = names(luad.exp)) 

# creates a new data.frame only with raw counts
mi <- luad.exp %>% select(-starts_with("xxx"))
# remove some unecessary characters of col.names
names(mi) = gsub(pattern = c(as.character("*.ead_count_")), 
                       replacement = "", x = names(mi)) 

# simplify col names even more
names(mi) = gsub(pattern = c(as.character("_01A.*|_01B.*|_01C.*")), replacement = "", x = names(mi)) #tirar tudo que tiver depois de 11A
names(mi)



# set the first column values back to row names so we can take means of duplicate columns
mi = tibble::column_to_rownames(mi, "miRNA_ID")


#### Take mean of duplicate columns, stay with these and delete the old ones
Adeno_miRNA <- as.data.frame( # sapply returns a list here, so we convert it to a data.frame
  sapply(unique(names(mi)), 
         function(col) rowMeans(mi[names(mi) == col]) 
  )
)


############### Sample stratifying ###############
setwd('D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/resultados/mirna_estratificacao')

# Stratify samples to get a 3:1 proportion between cancer and normal samples
vector <- colnames(Adeno_miRNA)
StrSample = list()
for (i in 1:10) {
  samples <- sample(vector, 150) # randomly choose the samples
  sample <- subset(Adeno_miRNA, select=c(samples)) # stratify it in new data.frames
  StrSample[[i]] <- sample # add it to your list
  names(StrSample)[[i]] <- paste0("Samples", i)
  dir.create(paste0("Samples",i, sep = ""))
  SampleNames = names(StrSample)
  setwd(file.path('D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/resultados/mirna_estratificacao/', SampleNames[[i]]))
  write.table(StrSample[[i]],
              paste0(SampleNames[i], ".txt", sep = ""), 
              sep = "\t",
              row.names = TRUE, 
              col.names = NA, 
              qmethod = c("escape", "double"))
  setwd('D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/resultados/mirna_estratificacao')
}
rm(i, sample, samples, vector)


############### Downloading TCGA Normal Lung Tissue miRNA expression files ############### 
queryNORMALm.exp <- GDCquery(project = "TCGA-LUAD", 
                            data.category = "Transcriptome Profiling",
                            data.type = "miRNA Expression Quantification",
                            workflow.type = "BCGSC miRNA Profiling",
                            experimental.strategy = "miRNA-Seq",
                            sample.type = c("Solid Tissue Normal"))
GDCdownload(queryNORMALm.exp)
NormalLung_m.exp <- GDCprepare(query = queryNORMALm.exp, save = TRUE, save.filename = "NormalLungExp.rda")


names(NormalLung_m.exp) = gsub(x = names(NormalLung_m.exp), pattern = "\\-", replacement = "_")
names(NormalLung_m.exp) = gsub(pattern = c(as.character("*.eads_per_million|*.ross_mapped")), 
                       replacement = "xxx", x = names(NormalLung_m.exp)) 

mi_n = NormalLung_m.exp %>% select(-starts_with("xxx"))

names(mi_n) = gsub(pattern = c(as.character("*.ead_count_")), 
                replacement = "", x = names(mi_n)) 

names(mi_n) = gsub(pattern = c(as.character("_11A.*|_11B.*|_11C.*")), replacement = "", x = names(mi_n)) #tirar tudo que tiver depois de 11A
names(mi_n)

mi_n = tibble::column_to_rownames(mi_n, "miRNA_ID")


# check for duplicates if you need too
Normal_miRNA <- as.data.frame( 
  sapply(unique(names(mi_n)), 
         function(col) rowMeans(mi_n[names(mi_n) == col]) 
  )
)

rm(queryLUAD.exp, queryNORMALm.exp, mi_n, mi)


############### save raw counts expression files ############### 
# if you want to manually explore it
write.table(Normal_miRNA, 
            file = 'normal.txt', 
            sep = "\t",
            row.names = TRUE, 
            col.names = NA, 
            qmethod = c("escape", "double"))

write.table(Adeno_miRNA, 
            file = 'Adeno.txt', 
            sep = "\t",
            row.names = TRUE, 
            col.names = NA, 
            qmethod = c("escape", "double"))


############### 