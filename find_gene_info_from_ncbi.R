##
## Finds the function/process/component from the gene ontology from NCBI
##

#Load libraries
require(rvest)

#Load original annotated datasets
fvb_conseq = read.csv("Output/fvb_conseq2.csv")
spret_conseq = read.csv("Output/spret_conseq2.csv")

##
## Load data from provean results and do general preprocessing
##
headers = c("#ROW_NO.", "INPUT","PROTEIN_ID","LENGTH", "STRAND", "CODON_CHANGE","POS","RESIDUE_REF","RESIDUE_ALT",
            "TYPE", "PRO_SCORE", "PREDICTION(cutoff=-2.5)","PRO_SEQ","#CLUSTER","SIFT_SCORE","PREDICTION(cutoff=0.05)",
            "MEDIAN_INFO","SIFT_SEQSEQ","dbSNP_ID")
fvb_sift = read.table("Input/fvb_v.result.tsv", sep="\t")
spret_sift = read.table("Input/spret_v.result.tsv", sep="\t")
names(fvb_sift) = headers
names(spret_sift) = headers

#Readd the chromosome and start pos back to the datasets for merging
add_identifiers = function(df){
  chr = regexpr("^[0-9]*", df$INPUT)
  pos = regexpr(",[0-9]*,", df$INPUT)
  
  #Recreate the seqnames column
  seqnames = substr(df$INPUT, chr, attr(chr, "match.length")) #Get the chromosome
  seqnames = paste0("chr", seqnames)
  
  #Recreate the start column
  start = substr(df$INPUT, pos + 1, attr(pos,"match.length") + 1) #Get start location, +1 for longer start
  start = sub(",", "", start)
  
  #Readd the columns back to the dataframe
  df$seqnames = seqnames
  df$start = start
  
  return(df)
}

#Parse out the chr and start site for merging with annotated data set
fvb_sift = add_identifiers(fvb_sift)
spret_sift = add_identifiers(spret_sift)

#Merge sift results with original annotated results
fvb = merge(fvb_conseq, fvb_sift, by=c("seqnames", "start"))
spret = merge(spret_conseq, spret_sift, by=c("seqnames", "start"))

##
## Create cleaned up dataset for eventual merging
##
fvb_data = fvb[,c("ID","GENENAME","seqnames", "GENEID", "TXID", "start", "end", "REF", "varAllele", "REFCODON", 
                  "VARCODON", "PROTEINLOC", "REFAA", "VARAA", "CONSEQUENCE","PREDICTION(cutoff=-2.5)", 
                  "PREDICTION(cutoff=0.05)")]
spret_data = spret[,c("ID","GENENAME","seqnames", "GENEID", "TXID", "start", "end", "REF", "varAllele", "REFCODON", "VARCODON",
                      "PROTEINLOC",  "REFAA", "VARAA", "CONSEQUENCE", "PREDICTION(cutoff=-2.5)", 
                      "PREDICTION(cutoff=0.05)")]

#Remove duplicates
fvb_data = unique(fvb_data)
spret_data = unique(spret_data)

#Convert ids with chr:loc into NA
fvb_data$ID[grep("^[0-9]*:", fvb_data$ID)] = NA
spret_data$ID[grep("^[0-9]*:", spret_data$ID)] = NA

#Probably add evidence column based on type vs prediction(whichever becomes delelerious)
fvb_data$Evidence = NA
spret_data$Evidence = NA

#If tolerated or neutral, non-damaging
fvb_data$Evidence[which(fvb_data$`PREDICTION(cutoff=-2.5)` == "Neutral" | fvb_data$`PREDICTION(cutoff=0.05)` == "Tolerated")] = "Non-Damaging"
spret_data$Evidence[which(spret_data$`PREDICTION(cutoff=-2.5)` == "Neutral" | spret_data$`PREDICTION(cutoff=0.05)` == "Tolerated")] = "Non-Damaging"

#If Synonymous and non-damaging -> synonymous
fvb_data$Evidence[which(fvb_data$CONSEQUENCE != "frameshift" & fvb_data$CONSEQUENCE != "nonsense" & 
                          fvb_data$CONSEQUENCE != "nonsynonymous" 
                        & fvb_data$Evidence == "Non-Damaging")] = "Synonymous"
spret_data$Evidence[which(spret_data$CONSEQUENCE != "frameshift" & spret_data$CONSEQUENCE != "nonsense" & 
                            spret_data$CONSEQUENCE != "nonsynonymous" & spret_data$Evidence == "Non-Damaging")] = "Synonymous"

#If a prediction is damaging, evidence = damaging
fvb_data$Evidence[which(fvb_data$`PREDICTION(cutoff=-2.5)` == "Deleterious" | fvb_data$`PREDICTION(cutoff=0.05)` == "Damaging")] = "Damaging"
spret_data$Evidence[which(spret_data$`PREDICTION(cutoff=-2.5)` == "Deleterious" | spret_data$`PREDICTION(cutoff=0.05)` == "Damaging")] = "Damaging"

#NA values have no evidence
fvb_data$Evidence[which(is.na(fvb_data$Evidence))] = "No Evidence"
spret_data$Evidence[which(is.na(spret_data$Evidence))] = "No Evidence"

#Frameshift/nonsense categorized as such
fvb_data$Evidence[which(fvb_data$CONSEQUENCE == "frameshift")] = "Frameshift"
spret_data$Evidence[which(spret_data$CONSEQUENCE == "frameshift")] = "Frameshift"
fvb_data$Evidence[which(fvb_data$CONSEQUENCE == "nonsense")] = "Nonsense"
spret_data$Evidence[which(spret_data$CONSEQUENCE == "nonsense")] = "Nonsense"

#Assign id to column based on evidence
fvb_data$rsIDs = NA
fvb_data$rsIDs =  NA

#Returns a character vector with ids that are considered damaging/frameshift/nonsense
find_damaging_ids = function(gene, dataset){
  gene_subset = subset(dataset, GENENAME == gene)
  id_list = (gene_subset$ID[gene_subset$Evidence == "Damaging" |
                              gene_subset$Evidence == "Frameshift" | 
                              gene_subset$Evidence == "Nonsense"])
  id_list = as.character(unique(id_list))
  
  #Checks for length of response and returns NA if nothing otherwise parses into single word if necessary
  if(length(id_list) != 0){
    if(length(id_list) > 1){
      return(paste(id_list,sep="", collapse="/"))
    } else{
      return(id_list)
    }
  } else{
    return("No Entries Found")
  }
}

#Finds and inputs the damaging ids into the column
for(gene in unique(fvb_data$GENENAME)){
  fvb_data$rsIDs[fvb_data$GENENAME == gene] = find_damaging_ids(gene, fvb_data)
}

for(gene in unique(spret_data$GENENAME)){
  spret_data$rsIDs[spret_data$GENENAME == gene] = find_damaging_ids(gene, spret_data)
}

#Save values of cleaned up dataset
write.csv(fvb_data, "Output/fvb_genes.csv", row.names = FALSE)
write.csv(spret_data, "Output/spret_genes.csv", row.names = FALSE)


##
## Create a dataframe between the two strains based on damage with possible comparsions
##
fvb_temp = fvb_data
spret_temp = spret_data
fvb_temp$Strain = "FVB"
spret_temp$Strain = "SPRET"
combined_strains = rbind(fvb_temp, spret_temp)

#Remove genes without any damaging aspects
combined_strains = combined_strains[!is.na(combined_strains$rsIDs),]

rm(fvb_temp, spret_temp)

#Order based on genename and start site
combined_strains = combined_strains[order(combined_strains$GENENAME, combined_strains$start),]

#Remove residual X
combined_strains$X = NULL

#Loops through and builds new data frame based on if there are duplicates for start/stop sites
single_var = combined_strains[1,]
for(i in 2:nrow(combined_strains)){
  f1 = single_var[nrow(single_var),]
  f2 = combined_strains[i,]
  f1$PROTEINLOC = NULL
  f2$PROTEINLOC = NULL
  testing = rbind(f1,f2)
  if(duplicated(testing)[2] == "FALSE"){ #Check if direct duplicates
    if(f1$start != f2$start && f1$end != f2$end){ #Check to see if start/stop the same place
        single_var = rbind(single_var, combined_strains[i,])        
    } else{
      if(f1$varAllele != f2$varAllele){ #Checks to see if those with same start or stop are different variant allele
        single_var = rbind(single_var, combined_strains[i,])        
      }
    }
  }
}

#Fix rownames for saving and reorder based on genename and start
combined_strains = combined_strains[order(combined_strains$GENENAME, combined_strains$start),]
rownames(combined_strains) = seq(1, nrow(combined_strains))

#Rename seqnames into chromosome and convert it into numeric
names(combined_strains)[which(names(combined_strains) == "seqnames")] = "Chrom"
combined_strains$Chrom = sub("chr", "", combined_strains$Chrom)
combined_strains$Chrom = as.numeric(combined_strains$Chrom)

#Convert to char to try to fix factor issue
combined_strains$rsIDs = as.character(combined_strains$rsIDs)

#Evidence IDs are based on gene not from gene/strain
for(gene in unique(combined_strains$GENENAME)){
  combined_strains$rsIDs[combined_strains$GENENAME == gene] = find_damaging_ids(gene, combined_strains)
}

#Create output summary dataframe
qtl_output = data.frame(
  GENENAME = unique(combined_strains$GENENAME)[1:length(unique(combined_strains$GENENAME)) - 1]
  , Chrom = NA, Evidence = NA ,rsIDs = NA
)

#Base url for querying ncbi
baseURL="http://www.ncbi.nlm.nih.gov/gene/"

#Add in to qtl_output, the highest level of evidence for the gene as well as the rsids
for (gene in unique(combined_strains$GENENAME)[1:length(unique(combined_strains$GENENAME)) - 1]) {
  gene_subset = subset(combined_strains, GENENAME == gene)
  
  qtl_output$Chrom[qtl_output$GENENAME %in% gene] = gene_subset$Chrom[1]
  
  qtl_output$rsIDs[qtl_output$GENENAME %in% gene] = gene_subset$rsIDs[1]
  if (length(gene_subset$Evidence[gene_subset$Evidence == "Nonsense"]) != 0) {
    qtl_output$Evidence[qtl_output$GENENAME %in% gene] = "Nonsense"
  }else if (length(gene_subset$Evidence[gene_subset$Evidence == "Frameshift"]) != 0) {
    qtl_output$Evidence[qtl_output$GENENAME %in% gene] = "Frameshift"
  }else if (length(gene_subset$Evidence[gene_subset$Evidence == "Damaging"]) != 0) {
    qtl_output$Evidence[qtl_output$GENENAME %in% gene] = "Damaging"
  }else if (length(gene_subset$Evidence[gene_subset$Evidence == "Non-damaging"]) != 0) {
    qtl_output$Evidence[qtl_output$GENENAME %in% gene] = "Non-damaging"
  } else{
    qtl_output$Evidence[qtl_output$GENENAME %in% gene] = "No Evidence"
  }
}

#Convert TXID into TXNAME
## Loads in the known genes for mouse
txdb_mm10 = TxDb.Mmusculus.UCSC.mm10.knownGene
txdb_mm10 = keepStandardChromosomes(txdb_mm10) #Remove weird chromosomes

## Get the TX Names
GRList = exonsBy(txdb_mm10)
tx_ids = names(GRList)
txnames = select(txdb_mm10, keys = tx_ids, columns = "TXNAME", keytype = "TXID")

## Merge with the txid
combined_strains = merge(txnames, combined_strains)

#Rearrange combined_strains
## ID GENENAME TXNAME Chrom Strain GENEID start end REF varAllele REFAA 
## VARAA PROTEINLOC CONSEQUENCE PREDICTION(cutoff..2.5) PREDICTION.cutoff(0.05) Evidence rsIDS
combined_strains = combined_strains[,c("ID", "GENENAME","TXNAME", "Chrom", "Strain", "GENEID", "start",
                                       "end", "REF", "varAllele", "REFAA", 
                                       "VARAA", "PROTEINLOC", "CONSEQUENCE", "PREDICTION(cutoff=-2.5)", 
                                       "PREDICTION(cutoff=0.05)", "Evidence", "rsIDs")]

#Output text files containing the summarized gene snp consequences as well as csv file containing original
write.table(qtl_output, "Output/annotation_summary.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.csv(combined_strains, "Output/annotation_data.csv", row.names = FALSE)

##
## Query NCBI for functional data for annotated genes
##

# Init. data for eventual saving
annot_functions = qtl_output
annot_functions$rsIDs = NULL
annot_functions$Functions = NA
annot_functions$Processes = NA
annot_functions$Components = NA

#Base url for querying ncbi
baseURL="http://www.ncbi.nlm.nih.gov/gene/"

# Find the functions of gene from ncbi
find_functions = function(gene_subset){
  gene.id = unique(gene_subset$GENEID)
  print(gene.id)
  
  gene_info = data.frame(Functions=NA, Processes=NA, Components=NA)
  if(is.na(gene.id)){
    return(NA)
  } else {
    link = paste0(baseURL, gene.id)
    link = read_html(link)
    
    #Parse out the function data, will take a while
    geneFunc <- link %>% html_nodes("td[headers = ont-func]") %>% html_text()
    geneFunc = gsub("\n", "", geneFunc)
    geneFunc = substring(geneFunc, 17)
    for(i in 1:length(geneFunc)){
      geneFunc[i] = substring(geneFunc[i], 1, nchar(geneFunc[i]) - 14)
    }
    
    #Parse out process data
    genePro <- link %>% html_nodes("td[headers = ont-comp]") %>% html_text()
    genePro = gsub("\n", "", genePro)
    genePro = substring(genePro, 17)
    for(i in 1:length(genePro)){
      genePro[i] = substring(genePro[i], 1, nchar(genePro[i]) - 14)
    }
    
    #Parse out component data
    geneComp <- link %>% html_nodes("td[headers = ont-proc]") %>% html_text()
    geneComp = gsub("\n", "", geneComp)
    geneComp = substring(geneComp, 17)
    for(i in 1:length(geneComp)){
      geneComp[i] = substring(geneComp[i], 1, nchar(geneComp[i]) - 14)
    }
    
    geneFunc = paste(geneFunc, collapse=", ")
    genePro = paste(genePro, collapse=", ")
    geneComp = paste(geneComp, collapse=", ")
    
    gene_info$Functions = geneFunc
    gene_info$Processes = genePro
    gene_info$Components = geneComp
    
    
    return(gene_info)
  }
}

# Subset out genes and find their function
for (gene in unique(combined_strains$GENENAME)[1:length(unique(combined_strains$GENENAME)) - 1]) {
  gene_subset = subset(combined_strains, GENENAME == gene)
  gene_info = find_functions(gene_subset)
  annot_functions$Functions[annot_functions$GENENAME == gene] = gene_info$Functions
  annot_functions$Processes[annot_functions$GENENAME == gene] = gene_info$Processes
  annot_functions$Components[annot_functions$GENENAME == gene] = gene_info$Components
}

annot_functions$Functions = gsub("_", " ", annot_functions$Functions)
annot_functions$Processes = gsub("_", " ", annot_functions$Processes)
annot_functions$Components = gsub("_", " ", annot_functions$Components)

write.csv(annot_functions, "Output/function_annot.csv", row.names = FALSE)

