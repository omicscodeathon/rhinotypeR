



PlotFrequency <- function(genotypeassigned,xlab='The assigned RV genotype',ylab='Sample genotype frequency'){
    # This function takes assignTypes() output and the user-supplied x and y-axis labels.
    # It tabulates the frequency of unique entries of the assigned_type column
    # Visualize using a bar plot the frequency of each genotype present in user samples
    plot=genotypeassigned %>%
    count(assigned_type) %>%
    ggplot(aes(x=assigned_type,y=n,fill=assigned_type))+
    geom_col()+
    labs(x=xlab, y=ylab)+
    theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5))
    return(plot)
}

# usage
PlotFrequency(data_type)



# Heatmap

library(gplots)
PlotPrototypeDistances <- function(pairwisedist, title="Dissimilarity Matrix Heatmap",xlabel='Prototype Refs',
                                  ylabel='Query sequences'){
    # This function takes a matrix of pairwise distance from allprotypeDistances func
    # Use gplot:heatmap.2, dark red shows more similarity
    # Modify the column names
    colnames(pairwisedist) <- sub(".*?_", "", colnames(pairwisedist)) # Comment out if you want to keep ref seq names as they are 
    heatmap <- heatmap.2(pairwisedist,
                         Rowv = NULL, Colv = NULL,  
                         dendrogram = 'none',       
                         trace = 'none',            
                         col = heat.colors(256),   
                         scale = 'row',             
                         margins = c(5, 10),        
                         main = title,
                        ylab=ylabel,
                        xlab=xlabel,
                        cex.axis = 0.8,
                        cex.lab = 0.8,
                        srtCol = 45,
                        adjCol = c(0, 0))   # Turn off this code to adjust the position
             
    return(heatmap$call) # Return the plot
}

# Usage
PlotPrototypeDistances(dist)