

# get input data from assignTypes.R

source("./scripts/04_assignTypes.R")



df <-assignTypes(pathToRef = "./data/output/RVBRefs.fasta", pathToQuery = "./data/testData.fasta", 
            "Tamura3p", 0.105)

  
  # Your provided list of values to randomly assign
replacement_values <- c("A1", "A21", "B13", "C11", "C53", "A101", "A57", "C3", "A87", "A9")

# Applying the random replacement
df$assigned_type <- ifelse(df$assigned_type == "unassigned",
                           sample(replacement_values, size = nrow(df), replace = TRUE),
                           df$assigned_type)

# Add 'species' column based on the first letter of 'assigned_type'
df$species <- substr(df$assigned_type, 1, 1)


# Load necessary library
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

# Aggregate data to summarize counts by species
species_counts <- aggregate(query ~ species, data = df, FUN = length)
species_counts$label <- paste("n =", species_counts$query)

# Plotting the circular bar plot
ggplot(species_counts, aes(x = species, y = query, fill = species)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(start = 0) + # Circular layout
  theme_void() +  # Remove background, gridlines, and text
  labs(title = "Species Distribution", x = NULL, y = NULL) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  theme(aspect.ratio = 1)


## 
# Type count

types_counts <- aggregate(query ~ assigned_type, data = df, FUN = length)
types_counts$label <- paste0(types_counts$assigned_type, ", ", types_counts$query)
types_counts$species <- substr(types_counts$assigned_type, 1, 1)
types_counts <- transform(types_counts, 
                          end_y = (query), 
                          start_y = c(0, head(cumsum(query), n=-1)))

ggplot(types_counts, aes(x = factor(assigned_type, levels = unique(assigned_type)), 
                         y = query, fill = species)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(start = 0) + # Circular layout
  theme_void() +  # Remove background, gridlines, and text
  labs(title = "Species Distribution", x = NULL, y = NULL) +
  geom_text(aes(y = end_y, label = label), vjust = -0.5, color = "black", size = 3) +
  theme(aspect.ratio = 1) # Ensure the plot is circular


# Simple bar chart
# Define colors for each species
color_map <- c("A" = "blue", "B" = "red", "C" = "green")
colors <- color_map[types_counts$species]

# Plot the barplot
barplot(height = types_counts$query, names.arg = types_counts$assigned_type, col = colors, 
        main = "Frequency of Types", xlab = "RV Type", ylab = "Count")

# Add legend
legend("topright", inset = c(0.01, -0.25), # Adjust the inset to move the legend
       legend = names(color_map), fill = color_map, title = "Species", 
       horiz = TRUE, bty = "n")

