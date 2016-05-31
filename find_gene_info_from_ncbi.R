##
## Finds the function/process/component from the gene ontology from NCBI
##

#Load libraries
require(rvest)

##
## Query NCBI for functional data for annotated genes
##

#Base url for querying ncbi
baseURL="http://www.ncbi.nlm.nih.gov/gene/"

# Find the functions of gene from ncbi
find_ontology = function(geneid){
    
    #Init return dataframe
    gene_info = data.frame(Functions=NA, Processes=NA, Components=NA)
    
    #Check for valid gene id
    if(is.na(geneid)){
        return(NA)
    } else {
        link = paste0(baseURL, geneid)
        link = read_html(link)
    
        #Parse out the function data, will take a while
        print("Parsing ontology - functions")
        geneFunc <- link %>% html_nodes("td[headers = ont-func]") %>% html_text()
        geneFunc = gsub("\n", "", geneFunc)
        geneFunc = substring(geneFunc, 17)
        for(i in 1:length(geneFunc)){
        geneFunc[i] = substring(geneFunc[i], 1, nchar(geneFunc[i]) - 14)
        }
    
        #Parse out process data
        print("Parsing ontology - processes")
        genePro <- link %>% html_nodes("td[headers = ont-comp]") %>% html_text()
        genePro = gsub("\n", "", genePro)
        genePro = substring(genePro, 17)
        for(i in 1:length(genePro)){
        genePro[i] = substring(genePro[i], 1, nchar(genePro[i]) - 14)
        }
    
        #Parse out component data
        print("Parsing ontology - components")
        geneComp <- link %>% html_nodes("td[headers = ont-proc]") %>% html_text()
        geneComp = gsub("\n", "", geneComp)
        geneComp = substring(geneComp, 17)
        for(i in 1:length(geneComp)){
        geneComp[i] = substring(geneComp[i], 1, nchar(geneComp[i]) - 14)
        }
    
        #Separate by comma
        geneFunc = paste(geneFunc, collapse=", ")
        genePro = paste(genePro, collapse=", ")
        geneComp = paste(geneComp, collapse=", ")
    
        #Assign values to return dataframe
        gene_info$Functions = geneFunc
        gene_info$Processes = genePro
        gene_info$Components = geneComp
    
        return(gene_info)
    }
}


