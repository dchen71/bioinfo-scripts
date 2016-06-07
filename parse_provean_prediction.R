##
## Function to parse Provean predictions
##

#Returns dataframe from parsed data containing anything in the family/domains section of uniprot
parse_provean_prediction = function(provean_file){
  #Load libraries
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)

  #Load data from provean results and do general preprocessing
  headers = c("#ROW_NO.", "INPUT","PROTEIN_ID","LENGTH", "STRAND", "CODON_CHANGE","POS","RESIDUE_REF","RESIDUE_ALT",
              "TYPE", "PRO_SCORE", "PREDICTION(cutoff=-2.5)","PRO_SEQ","#CLUSTER","SIFT_SCORE","PREDICTION(cutoff=0.05)",
              "MEDIAN_INFO","SIFT_SEQSEQ","dbSNP_ID")
  provean_read = read.table(provean_file, sep="\t")
  names(provean_read) = headers
  
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
  provean_read = add_identifiers(provean_read)
  
  return(provean_read)
}