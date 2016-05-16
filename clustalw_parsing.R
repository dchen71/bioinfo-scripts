##
## Parse ClustalW results and return back dataframe with base amd -1, 0, 1, 2 based on conservation
##

#Read input
raw_data = read.table("muscle-I20160516-190153-0407-29699783-pg.clw", sep="\t")

#Parses clustalW results and return back dataframe
parse_clustal = function(clustal_results, num_fasta = 2){
    #Return results for all entries
    get_results = function(entry){
        results = clustal_results[seq(entry + 1,nrow(clustal_results), by= num_fasta + 1),]
        results = substr(results, 27, nchar(as.character(results[1])))
        results = paste0(results, sep="", collapse="") #Single line
    }
    
    #Get back single line version of all entries by metaprogramming based on number of fasta inputs
    get_reads <- function() {
        for(i in 1:num_fasta + 1) {
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
    con = get(paste0("result.",num_fasta + 1))
    
    #Loop to append graph df
    for(i in 1:nchar(con)){
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
        graph_df$MUS[i] = substr(mus,i,i)
        graph_df$HUMAN[i] = substr(human,i,i)
    }
    
    #return(graph_df)
}




get_reads(5)

genMyFuns(c(7, 9, 10))


#Plot results
write.csv(graph_df, "../Input/parsed_con.csv", row.names = FALSE)
