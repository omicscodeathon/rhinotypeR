# input data from assignTypes.R


plotFrequency <- function(assignedTypesDF, showLegend=FALSE){
  # Add 'species' column based on the first letter of 'assigned_type'
  assignedTypesDF$species <- substr(assignedTypesDF$assignedType, 1, 1)
  
  
  types_counts <- aggregate(query ~ assignedType, data = assignedTypesDF, FUN = length)
  types_counts$label <- paste0(types_counts$assignedType, ", ", types_counts$query)
  types_counts$species <- substr(types_counts$assignedType, 1, 1)
  types_counts <- transform(types_counts,
                            end_y = (query),
                            start_y = c(0, head(cumsum(query), n=-1)))
  
  types_counts$species[types_counts$species == "u"] <- "Other"
  
  # Simple bar chart
  # Define colors for each species
  color_map <- c("A" = "blue", "B" = "red", "C" = "green", "Other" = "grey")
  colors <- color_map[types_counts$species]
  
  # Plot the barplot
  barplot(height = types_counts$query, names.arg = types_counts$assignedType, col = colors, 
          main = "Frequency of Types", xlab = "RV Type", ylab = "Count")
  if (showLegend){
    # Add legend
    legend("topright", inset = c(0.01, -0.17), # Adjust the inset to move the legend
           legend = names(color_map), fill = color_map, title = "Species", 
           horiz = TRUE, bty = "n")
  }
  
}
