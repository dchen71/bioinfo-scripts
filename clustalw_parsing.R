##
## Parse ClustalW results and return back dataframe with base amd -1, 0, 1, 2 based on conservation
##

#Parses clustalW results and return back dataframe
parse_clustal = function(clustal_results, num_fasta = 2){
    #Return results for all entries
    get_results = function(entry){
        results = clustal_results[seq(entry + 1,nrow(clustal_results), by= num_fasta + 1),]
        results = substr(results, 27, nchar(as.character(results[1])))
        results = paste0(results, sep="", collapse="") #Single line
        return(results)
    }
    
    #Get back single line version of all entries by metaprogramming based on number of fasta inputs
    get_reads <- function() {
        for(i in 1:(num_fasta + 1)) {
            fName <- paste("results.", i, sep="")
            assign(fName, eval(
                substitute(
                    get_results(i)
                )
            ),
            envir=parent.frame()
            )
        }
    }
    
    get_reads()
    
    #Create data frame for graphing
    graph_df = data.frame(POS = seq(1, nchar(results.1)), CON = NA)
    
    #Init empty rows based on number of fasta inputs
    #New rows based on first fasta file to last
    graph_df = cbind(graph_df, data.frame(matrix(ncol=num_fasta)))
    
    #Get the clustal results
    con = get(paste0("results.",num_fasta + 1))
    
    #Loop to append graph df
    for(i in 1:nchar(con)){
        #Code based on clustal grading
        # * = 0
        # : = 2
        # . = 1
        # space = -1
        curr_char = substr(con,i,i)
        if(curr_char == "*"){
            graph_df$CON[i] = 0
        } else if(curr_char == ":"){
            graph_df$CON[i] = 2
        } else if(curr_char == "."){
            graph_df$CON[i] = 1
        } else{
            graph_df$CON[i] = -1
        }
        #Adds fasta character
        for(fasta_num in 1:num_fasta){
            #Loop through to move through dataframe columns
            graph_df[i,(2 + fasta_num)] = substr(get(paste0("results.",fasta_num)),i,i)

        }

    }
    
    return(graph_df)
}