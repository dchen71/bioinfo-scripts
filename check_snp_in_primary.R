##
## Checks if a given rsid are included in primary transcript or not from dbSNP
## Requires: rsid, gene id, protein id, mrna id, and contig id from dbSNP
##

# Find if the rsid is included in the main transcript of gene from dbsnp
find_canonical = function(rsid, geneid, protein_id, mrna_id, contig_id){
  
  # Updating link for the normal transcript
  link = "http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?locusId="
  link = paste0(link, geneid)
  link = paste0(link, "&ctg=", contig_id)
  link = paste0(link, "&mrna=", mrna_id)
  link = paste0(link, "&prot=", protein_id)
  link = read_html(link)
  
  #Collects rsids
  rsids <- link %>% html_nodes(".in_gm_rs") %>% html_text()
  
  if(!is.na(grep(rsid, rsids)[1])){
    return("Primary")
  }
  else{
    return("Variant")
  }
}
