##
## Function to find domain data per gene UniprotID
##

#Returns dataframe from parsed data containing anything in the family/domains section of uniprot
find_uniprot_domains = function(uniprot_id){
  #Load libraries
  require(rvest)

  #Base url for querying uniprot
  baseURL="http://www.uniprot.org/uniprot/"
  
  # Link to mtss1 main
  link = paste0(baseURL, uniprot_id)
  link = read_html(link)
  
  #Collect the gene name
  gene_name <- link %>% html_nodes("#content-gene h2") %>% html_text()
  
  #Collect the species
  organism_name <- link %>% html_nodes("#content-organism em") %>% html_text()
  
  #Collect domain info
  domain_data <- link %>% html_nodes("#family_and_domains .featureTable .feature_row td") %>% html_text()
  
  #Check if no domain
  if(length(domain_data) == 0){ 
    return(NA)
  } else {
    features = domain_data[seq(1,length(domain_data), 7)]
    positions = domain_data[seq(2,length(domain_data), 7)]
    descriptions = domain_data[seq(4,length(domain_data), 7)]
    
    #Clean up features
    for(i in 1:length(features)){
      features[i] = substr(features[i], 1, nchar(features[i]) - 1)
      features[i] = substr(features[i], gregexpr("</p>", features[i])[[1]][2] + 4, nchar(features[i]))
    }
    
    #Clean up descriptions
    for(i in 1:length(descriptions)){
      prosite_index = as.numeric(regexpr("PROSITE", descriptions[i]))
      interpro_index = as.numeric(regexpr("InterPro", descriptions[i]))
      seq_index = as.numeric(regexpr("Sequence analysis", descriptions[i]))
      similarity_index = as.numeric(regexpr("By similarity", descriptions[i]))
      if(prosite_index != -1){ #If prosite aint there
        descriptions[i] = substr(descriptions[i], 1, prosite_index - 1)
      } else if(interpro_index != -1){
        descriptions[i] = substr(descriptions[i], 1, interpro_index - 1)
      } else if(similarity_index != -1){
        descriptions[i] = substr(descriptions[i], 1, similarity_index - 1)
      } else if(seq_index != -1){
        descriptions[i] = NA
      }
    }
    
    #Lazy init char vector for start and end positions
    start = positions
    end = positions
    
    #Split positions into start and end
    for(i in 1:length(positions)){
      start_length = attr(regexpr("[0-9]*", positions[i]), "match.length")
      start[i] = substr(positions[i], 1, start_length)
      
      end[i] = substr(positions[i], start_length + 4, nchar(positions[i]))
    }
    
    #Collect Uniprot protein status
    reviewed <- link %>% html_nodes("#content-status a") %>% html_text()
    
    #Find how uniprot evaluated protein
    evidence <- link %>% html_nodes("#content-status .context-help") %>% html_text()
    evidence = evidence[2]
    evidence = substr(evidence, 27,as.numeric(regexpr(" level", evidence) - 1))
    if(!is.na(evidence) & nchar(evidence) == 0){
      evidence = "Inferred"
    }
    
    #Create data frame to return
    transcript_data = data.frame(GENENAME = gene_name,
                                 UniprotID = uniprot_id, 
                                 Organism = organism_name,
                                 Reviewed = reviewed, 
                                 Evidence = evidence, 
                                 Features = features, 
                                 Start = start,
                                 End = end,
                                 Descriptions = descriptions)
    
    for(i in 1:nrow(transcript_data)){
      if(is.na(transcript_data$Descriptions[i])){
        next
      } else{
        similarity_index = as.numeric(regexpr("By similarity", transcript_data$Descriptions[i]))
        if(similarity_index != -1){
          transcript_data$Descriptions[i] = (substr(transcript_data$Descriptions[i], 1, similarity_index - 1))
        }
      }
    }
    
    return(transcript_data)
  }
}