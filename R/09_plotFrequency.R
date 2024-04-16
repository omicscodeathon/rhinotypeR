source("R/04_assignTypes.R")
source("R/02_readFasta.R")

#target_fasta_2 <- "./data/testData.fasta"
#usethis::use_data(target_fasta_2)
#refSeq <- "./data/output/RVBRefs.fasta"
#usethis::use_data(refSeq)

plotFrequency <- function(target_data, ref_data, model = "Tamura3p") {
  
  df <- assignTypes(ref_data, readFasta(target_data), "Tamura3p", 0.105)
  
    # Your provided list of values to randomly assign
  replacement_values <- c("A1", "A21", "B13", "C11", "C53", "A101", "A57", "C3", "A87", "A9")
  
  # Applying the random replacement
  df$assigned_type <- ifelse(df$assigned_type == "unassigned",
                             sample(replacement_values, size = nrow(df), replace = TRUE),
                             df$assigned_type)
  
  # Add 'species' column based on the first letter of 'assigned_type'
  df$species <- substr(df$assigned_type, 1, 1)
  
  
  types_counts <- aggregate(query ~ assigned_type, data = df, FUN = length)
  types_counts$label <- paste0(types_counts$assigned_type, ", ", types_counts$query)
  types_counts$species <- substr(types_counts$assigned_type, 1, 1)
  types_counts <- transform(types_counts,
                            end_y = (query),
                            start_y = c(0, head(cumsum(query), n=-1)))

  # Simple bar chart
  # Define colors for each species
  color_map <- c("A" = "blue", "B" = "red", "C" = "green")
  colors <- color_map[types_counts$species]
  
  # Plot the barplot
  barplot(height = types_counts$query, names.arg = types_counts$assigned_type, col = colors, 
          main = "Frequency of Types", xlab = "RV Type", ylab = "Count")
  
  # Add legend
  legend("topright", inset = c(0.01, -0.17), # Adjust the inset to move the legend
         legend = names(color_map), fill = color_map, title = "Species", 
         horiz = TRUE, bty = "n")

}
