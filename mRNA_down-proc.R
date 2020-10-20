
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


############### Download TCGA Lung cancer mRNA expression files ###############

setwd('D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/resultados/mrna_estratificacao') # set you working directory

TCGAbiolinks:::getGDCprojects()$project_id  # gives a list of GDC projects IDs according to cancer/sample/data source types
TCGAbiolinks:::getProjectSummary("TCGA-LUAD") # gives a list with possible data categories of files of lung adenocarcinoma (LUAD) patients from TCGA that can be downloaded

# first you set the query arguments to download mRNA counts files of LUAD samples
queryLUAD.exp <- GDCquery(project = "TCGA-LUAD", # choose a TCGA project of interest
                          data.category = "Transcriptome Profiling", # (Biospecimen, Clinical, CNV, DNA methylation, Sequencing Reads, SNV, transcriptome profiling, Gene expression, Protein expression, Raw microarray data, Raw sequencing data)
                          data.type = "Gene Expression Quantification", # data type to filter the files to download (Gene Expression Quantification, Isoform Expression Quantification, miRNA Expression Quantification, Splice Junction Quantification)
                          workflow.type = "HTSeq - Counts", # specify what type of counts you want (HTSeq-Counts, HTSeq-FPKM, HTSeq-FPKM-UQ, STAR-Counts)
                          experimental.strategy = "RNA-Seq", #Filter to experimental strategy (WXS, RNA-Seq, miRNA-Seq, Genotyping Array, DNA-Seq, Methylation array, Protein expression array, CGH array, VALIDATION, Gene expression array, WGS, MSI-Mono-Dinucleotide Assay, miRNA expression array, Mixed strategies, AMPLICON, Exon array, Total RNA-Seq, Capillary sequencing, Bisulfite-Seq)
                          sample.type = c("Primary Tumor")) # you can find sample type complete list on https://www.bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/query.html#sampletype_options

GDCdownload(queryLUAD.exp) # then download the data to your computer (can be found in the active working directory)
luad.exp <- GDCprepare(query = queryLUAD.exp, save = TRUE, save.filename = "luadExp.rda") # and import it to your R enviroment
# The data downloaded here is from TCGA database. The LUAD project has data of lung adenocarcinoma tissue samples. We downloaded gene expression RNA-seq data samples from 533 patients with (Primary solid Tumor) lung adenocarcinoma.


# Remove ffpe samples if you need/want to (here it removes 11)
sem_ffpe <- luad.exp[, !luad.exp$is_ffpe]

# Filter low correlated samples
AD_UnprocBarc <- TCGAanalyze_Preprocessing(object = sem_ffpe,
                                           cor.cut = 0.8, # choose the correlation cut-off
                                           datatype = "HTSeq - Counts", #specify the counts type
                                           filename = "luad_filt_ffpe_IlluminaHiSeq_RNASeqV2.png") # choose .png file name for the correlation heatmap output
# This will also extract a matrix only with the raw counts from previous downloaded data

# change matrix to dataframe so we can use gsub to modify columns names
AD_UnprocBarc = as.data.frame(AD_UnprocBarc)

# Adjust columns names exchanging hyphens to underlines and remove some unecessary characters col.names
names(AD_UnprocBarc) = gsub(x = names(AD_UnprocBarc), pattern = "\\-", replacement = "_")
names(AD_UnprocBarc) = gsub(pattern = c(as.character("_01A.*|_01B.*|_01C.*")), 
                           replacement = "", x = names(AD_UnprocBarc)) 

# Take mean of duplicate columns, stay with these and delete the old ones
Adeno_RNA <- as.data.frame( # sapply returns a list here, so we convert it to a data.frame
  sapply(unique(names(AD_UnprocBarc)), 
         function(col) rowMeans(AD_UnprocBarc[names(AD_UnprocBarc) == col]) 
  )
)



############### Sample stratifying ###############
setwd('D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/resultados/mrna_estratificacao')

# Stratify samples to get a 3:1 proportion between cancer and normal samples
vector <- colnames(Adeno_RNA)
StrSample = list()
for (i in 1:10) {
  samples <- sample(vector, 150) # randomly choose the samples
  sample <- subset(Adeno_RNA, select=c(samples)) # stratify it in new data.frames
  StrSample[[i]] <- sample # add it to your list
  names(StrSample)[[i]] <- paste0("Samples", i)
  dir.create(paste0("Samples",i, sep = ""))
  SampleNames = names(StrSample)
  setwd(file.path('D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/resultados/mrna_estratificacao/', SampleNames[[i]]))
  write.table(StrSample[[i]],
              paste0(SampleNames[i], ".txt", sep = ""), 
              sep = "\t",
              row.names = TRUE, 
              col.names = NA, 
              qmethod = c("escape", "double"))
  setwd('D:/Faculdade/pibit/PIBIT/Finalização - pibit/merge/resultados/mrna_estratificacao')
}
rm(i, sample, samples, vector)



############### Download TCGA Normal Lung Tissue mRNA expression files ###############

# First you set the query arguments to download mRNA counts files of normal tissue samples
queryNORMAL.exp <- GDCquery(project = "TCGA-LUAD", 
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification",
                            workflow.type = "HTSeq - Counts",
                            experimental.strategy = "RNA-Seq",
                            sample.type = c("Solid Tissue Normal"))
GDCdownload(queryNORMAL.exp)
NormalLung.exp <- GDCprepare(query = queryNORMAL.exp, save = FALSE, save.filename = "NormalLungExp.rda")
# The data downloaded here is from TCGA database. The LUAD project has data of normal lung tissue samples too. We downloaded gene expression RNA-seq data samples from patients with (Primary solid Tumor) lung adenocarcinoma.

# This time we don't have ffpe samples to filter out but you can do it if you need

# Now repeat some of the steps we did before
NO_UnprocBarc <- TCGAanalyze_Preprocessing(object = NormalLung.exp,
                                                 cor.cut = 0.8,    
                                                 datatype = "HTSeq - Counts",
                                                 filename = "normal_filt_ffpe_IlluminaHiSeq_RNASeqV2.png")

# change matrix to dataframe so we can use gsub to modify columns names
NO_UnprocBarc = as.data.frame(NO_UnprocBarc)

# Adjust columns names exchanging hyphens to underlines and remove some unecessary characters col.names
names(NO_UnprocBarc) = gsub(x = names(NO_UnprocBarc), pattern = "\\-", replacement = "_")
names(NO_UnprocBarc) = gsub(pattern = c(as.character("_11A.*|_11B.*|_11C.*")), # here the pattern to be substituted is different
                         replacement = "", x = names(NO_UnprocBarc)) 
names(NO_UnprocBarc) # check

# Take mean of duplicate columns, stay with these and delete the old ones
Normal_RNA <- as.data.frame( # sapply returns a list here, so we convert it to a data.frame
  sapply(unique(names(NO_UnprocBarc)), 
         function(col) rowMeans(NO_UnprocBarc[names(NO_UnprocBarc) == col]) 
  )
) # take mean if you need too... but this time we don't have duplicates



############### Remove the excess of files before continue with analysis
rm(sem_ffpe, queryLUAD.exp, queryNORMAL.exp, NO_UnprocBarc, AD_UnprocBarc)
