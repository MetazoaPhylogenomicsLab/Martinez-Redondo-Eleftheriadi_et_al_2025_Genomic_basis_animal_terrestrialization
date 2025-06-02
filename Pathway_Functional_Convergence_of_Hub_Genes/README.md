# Pathways Functional Convergence

---

To test functional convergence at the pathway level we assigned KEGG pathways using eggNOG-mapper v2.1.12 to each hub gene involved in the respective stress responses. A pathway was considered relevant to the response if it comprised at least 1% of the total hub genes for that specific species-stress context.

We have summarised the results of eggNOG-mapper in this file summaryAllSpAllExp.RData

which can be found in the figshare Pathway_Functional_Convergence_of_Hub_Genes directory

First, we the pathways that are present in more then 75% of the species, combining all experiments:

[pathways_plots.R](Pathways%20Functional%20Convergence%20203f4c09425a803f946dca91f0da54dc/pathways_plots.r)

```r
# Set working directory
setwd("./Pathway_analysis")

# Load required libraries
library(dplyr)
library(stringr)
library(tibble)
library(pbmcapply)
library(Biostrings)
library(pheatmap)
library(KEGGREST)

################################################################################
# Step 1: Build KEGG pathway hierarchy
################################################################################

brite_raw <- keggGet("br:br08901")[[1]]
lines <- str_split(brite_raw, "\n")[[1]]

records <- list()
level1 <- level2 <- NULL

for (line in lines) {
  line <- str_trim(line)
  if (str_starts(line, "A")) {
    level1 <- str_remove(line, "^A\\s*")
  } else if (str_starts(line, "B")) {
    level2 <- str_remove(line, "^B\\s*")
  } else if (str_starts(line, "C")) {
    path_id <- str_extract(line, "\\d{5}")
    path_id <- paste0("path:map:", path_id)
    path_name <- str_trim(str_remove(line, "^C\\s*\\d{5}\\s*"))
    if (!is.na(path_id)) {
      records[[length(records) + 1]] <- tibble(
        Pathway_ID = path_id,
        Level3 = path_name,
        Level2 = level2,
        Level1 = level1
      )
    }
  }
}

# Combine and filter out non-relevant categories
pathway_df <- bind_rows(records)
pathway_df <- pathway_df[!pathway_df$Level1 %in% c("Human Diseases", "Drug Development"), ]

################################################################################
# Step 2: Load summarized annotation data
################################################################################

load("summaryAllSpAllExp.RData")

################################################################################
# Step 3: Build presence/absence matrices per experiment
################################################################################

# Pathway list
all_pathways <- sub("path:map:", "map", pathway_df$Pathway_ID)

# Keep all desired species (without collapsing TLON/TMEL)
desired_column_order <- c("CELE", "LMAR", "PEAU", "PLAE", "MISO", "ONUN", "SMED", "PTUR",
                          "SPEC", "TPIS", "PACU", "COCO", "EAND", "HMED", "LEPN", "TLON", "TMEL")

# List of experiments
exp_list <- c("CD", "CF", "Hyper", "Hypo", "OR", "UV0", "UV15D", "UV15L", "UV24D", "UV24L", "VL")

# Initialize list to hold all matrices
heatMaps_byExp_full <- list()

# Loop through experiments
for (exp_name in exp_list) {
  message("Processing: ", exp_name)
  
  countsbyExp <- lapply(summaryAllSpAllExp, function(df) {
    if (exp_name %in% colnames(df)) {
      relevant_rows <- df[(!df$pathways %in% c("", "-")) & df[[exp_name]], ]
      sapply(all_pathways, function(path_id) {
        sum(grepl(path_id, relevant_rows$pathways)) / sum(df[[exp_name]]) > 0.01
      })
    } else {
      rep(FALSE, length(all_pathways))
    }
  })
  
  # Combine into matrix
  countsbyExp <- do.call(cbind, countsbyExp)
  colnames(countsbyExp) <- sub("SPEC1", "SPEC", names(summaryAllSpAllExp))
  
  # Reorder columns to match desired order
  common_species <- intersect(desired_column_order, colnames(countsbyExp))
  countsbyExp <- countsbyExp[, common_species, drop = FALSE]
  rownames(countsbyExp) <- all_pathways
  
  # Store in output list
  heatMaps_byExp_full[[exp_name]] <- countsbyExp
}

################################################################################
# Step 4: Save matrices (optional)
################################################################################

# Save each matrix as CSV (presence/absence as 1/0)
for (exp in names(heatMaps_byExp_full)) {
  write.csv(heatMaps_byExp_full[[exp]] * 1, paste0("PresenceAbsence_", exp, ".csv"), row.names = TRUE)
}

################################################################################
#############    Make a plot by combining all experiments      #################
################################################################################

# Step 1: Combine all experiment matrices into one big matrix (species_experiment columns)
# Combine all experiment matrices into one big matrix (species_experiment columns)
combinedMatrix <- do.call(cbind, lapply(names(heatMaps_byExp_full), function(exp_name) {
  mat <- heatMaps_byExp_full[[exp_name]]
  colnames(mat) <- paste(colnames(mat), exp_name, sep = "_")
  return(mat)
}))

# Extract species names
all_colnames <- colnames(combinedMatrix)
base_species <- unique(sub("_.*", "", all_colnames))

# Replace TLON and TMEL with TETR
base_species <- unique(sub("TLON|TMEL", "TETR", base_species))

# Merge experiments and species (including TLON + TMEL → TETR)
mergedMatrix <- sapply(base_species, function(sp) {
  if (sp == "TETR") {
    # Combine TLON and TMEL
    cols <- grep("^(TLON|TMEL)_", all_colnames)
  } else {
    cols <- grep(paste0("^", sp, "_"), all_colnames)
  }
  rowSums(combinedMatrix[, cols, drop = FALSE]) > 0  # TRUE if any experiment is TRUE
})

mergedMatrix_binary <- mergedMatrix * 1

# Save to CSV
write.csv(mergedMatrix_binary, "PresenceAbsence_MergedBySpecies_TETR.csv", row.names = TRUE)
# Step 6: Keep only pathways present in >75% of species
threshold <- 0.75
species_count <- ncol(mergedMatrix)
pathway_freq <- rowSums(mergedMatrix) / species_count

# Filter pathways
filteredMatrix <- mergedMatrix[pathway_freq > threshold, ]

# Optional: convert to binary
filteredMatrix_binary <- filteredMatrix * 1

# # Match full pathway names using pathway_df
# pathway_ids <- rownames(filteredMatrix)
# full_names <- sapply(pathway_ids, function(pid) {
#   match_idx <- which(sub("path:map:", "map", pathway_df$Pathway_ID) == pid)
#   if (length(match_idx) > 0) {
#     pathway_df$Level3[match_idx]
#   } else {
#     pid  # fallback: use ID if name not found
#   }
# })
# 
# # Update rownames to full pathway names
# rownames(filteredMatrix) <- full_names
# filteredMatrix_binary <- filteredMatrix * 1

# Save filtered matrix
#write.csv(filteredMatrix_binary, "PresenceAbsence_75percSpecies_TETR.csv", row.names = TRUE)

# Step 1: Extract pathway IDs from the filtered matrix 
filtered_pathway_ids <- rownames(filteredMatrix)

# Step 2: Find corresponding pathway IDs in the original per-experiment matrices
# Create a mapping from full names to pathway IDs (this assumes you already have `pathway_df`)
pathway_id_map <- setNames(pathway_df$Pathway_ID, pathway_df$Level3)

# Step 3: For each experiment, filter pathways using the mapped IDs 
# So each experiment has the samepathways with the merged heatmap
heatMaps_byExp_filtered <- list()
heatMaps_byExp_full
for (exp_name in names(heatMaps_byExp_full)) {
  # Get the matrix for the current experiment
  exp_matrix <- heatMaps_byExp_full[[exp_name]]
  
  # Map full names to pathway IDs
  pathway_ids_in_exp <- rownames(exp_matrix)
                               
  # Filter the pathways by the IDs present in filteredMatrix_binary
  # First, get the subset of pathway IDs that are in both the experiment matrix and the filtered matrix
  matching_pathway_ids <- intersect(filtered_pathway_ids, pathway_ids_in_exp)
  
  # Subset the experiment matrix by matching pathway IDs
  exp_matrix_filtered <- exp_matrix[matching_pathway_ids, , drop = FALSE]
  
  # Store the filtered matrix in the list
  heatMaps_byExp_filtered[[exp_name]] <- exp_matrix_filtered
}
  
# Now `heatMaps_byExp_filtered` contains the filtered matrices with only the pathways present in `filteredMatrix_binary`

# plot heatmaps for each experiment in heatMaps_byExp_filtered
for (exp_name in names(heatMaps_byExp_filtered)) {
  # Get the filtered matrix for the current experiment
  exp_matrix_filtered <- heatMaps_byExp_filtered[[exp_name]]

  # Plot the heatmap ******
  pheatmap(exp_matrix_filtered * 1,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,  # Show the pathway names
           show_colnames = TRUE,  # Show the species/experiment names
           color = c("#ECEAE6", "#005f73"),  # Use binary colors (white for 0, black for 1)
           main = paste("Pathway Presence -", exp_name),
           #border_color = "grey",  # Optional: Add border for clarity
           cellwidth = 10,  # Optional: Adjust cell width
           cellheight = 10,  # Optional: Adjust cell height
           fontsize = 8
           )  # Optional: Adjust font size for row/column names
  dev.off()

}
# Save the filtered matrices to CSV for each experiment
#for (exp_name in names(heatMaps_byExp_filtered)) {
 # write.csv(heatMaps_byExp_filtered[[exp_name]], paste0("Filtered_PresenceAbsence_", exp_name, ".csv"), row.names = TRUE)
#}

# # Match full pathway names using pathway_df
pathway_ids <- rownames(filteredMatrix)
full_names <- sapply(pathway_ids, function(pid) {
  match_idx <- which(sub("path:map:", "map", pathway_df$Pathway_ID) == pid)
  if (length(match_idx) > 0) {
    pathway_df$Level3[match_idx]
  } else {
    pid  # fallback: use ID if name not found
  }
})

# Update rownames to full pathway names
rownames(filteredMatrix) <- full_names
filteredMatrix_binary <- filteredMatrix * 1
pdf("heatmap_allExperiments.pdf", width = 10, height = 10)  # Adjust size as needed
pheatmap(filteredMatrix_binary,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("#ECEAE6", "#005f73"))(2),
         cellwidth = 10,  # Optional: Adjust cell width
         cellheight = 10,  # Optional: Adjust cell height
         fontsize = 8,
         main = "Pathways in >75% of Species")
dev.off()

# Define a function to map pathway IDs to full pathway names using pathway_df
map_full_names <- function(pathway_ids) {
  full_names <- sapply(pathway_ids, function(pid) {
    match_idx <- which(sub("path:map:", "map", pathway_df$Pathway_ID) == pid)
    if (length(match_idx) > 0) {
      pathway_df$Level3[match_idx]
    } else {
      pid  # fallback: use ID if name not found
    }
  })
  return(full_names)
}

# Loop through each experiment's filtered matrix and update row names. This will generate the full set of heatmaps, however for the main figure only keep CD, Hyper, OR, UV15L and VL
first_plot <- TRUE  # Flag to control row name display

for (exp_name in names(heatMaps_byExp_filtered)) {
  # Get the filtered matrix for the current experiment
  exp_matrix_filtered <- heatMaps_byExp_filtered[[exp_name]]
  
  # Get the pathway IDs (rownames) in the experiment matrix
  pathway_ids <- rownames(exp_matrix_filtered)
  
  # Map these IDs to full pathway names using the map_full_names function
  full_names <- map_full_names(pathway_ids)
  
  # Update the rownames of the experiment matrix with the full pathway names
  rownames(exp_matrix_filtered) <- full_names
  
  # Convert the matrix to binary (0/1)
  exp_matrix_filtered_binary <- exp_matrix_filtered * 1
  
  # Remove columns with species without data in a given experiment from the matrix
  keep <- sapply(colnames(exp_matrix_filtered_binary), function(x){
    exp_name %in% names(summaryAllSpAllExp[[x]])
  })
  exp_matrix_filtered_binary <- exp_matrix_filtered_binary[, which(keep)]
  
  # Join the TLON and TMEL columns
  if(sum(colnames(exp_matrix_filtered_binary) %in% c("TLON", "TMEL")) > 0){
    if(sum(colnames(exp_matrix_filtered_binary) %in% c("TLON", "TMEL")) > 1){
      newVec <- rowSums(exp_matrix_filtered_binary[, which(colnames(exp_matrix_filtered_binary) %in% c("TLON", "TMEL"))]) > 0
      exp_matrix_filtered_binary <- cbind(exp_matrix_filtered_binary, newVec)
      colnames(exp_matrix_filtered_binary)[ncol(exp_matrix_filtered_binary)] <- "TETR"
      exp_matrix_filtered_binary <- exp_matrix_filtered_binary[, which(!colnames(exp_matrix_filtered_binary) %in% c("TLON", "TMEL"))]
    } else{
      colnames(exp_matrix_filtered_binary)[which(colnames(exp_matrix_filtered_binary) %in% c("TLON", "TMEL"))] <- "TETR"
    }
  }
  # Define PDF filename
  pdf_filename <- paste0("heatmap_", exp_name, ".pdf")
  # Create and save the heatmap
  pdf(pdf_filename, width = 8, height = 8)  # Adjust size as needed
  pheatmap(exp_matrix_filtered_binary,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = first_plot,  # Only TRUE for the first plot
           show_colnames = TRUE,
           color = c("#ECEAE6", "#005f73"),
           main = paste("Pathway Presence -", exp_name),
           legend=FALSE,
           #border_color = "grey",
           cellwidth = 8,
           cellheight = 8,
           fontsize = 8)
  dev.off()
  
  # After the first plot, turn off rownames
  first_plot <- FALSE
}
```

Let’s plot now the presence/absence heatmap of each experiment without filtering the pathways as before:

```r
# Loop through each experiment's full matrix and update row names. This will generate the full set of heatmaps
for (exp_name in names(heatMaps_byExp_full)) {
  # Get the matrix for the current experiment
  exp_matrix <- heatMaps_byExp_full[[exp_name]]
  
  # Get the pathway IDs (rownames) in the experiment matrix
  pathway_ids <- rownames(exp_matrix)
  
  # Map these IDs to full pathway names using the map_full_names function
  full_names <- map_full_names(pathway_ids)
  
  # Update the rownames of the experiment matrix with the full pathway names
  rownames(exp_matrix) <- full_names
  
  # Convert the matrix to binary (0/1)
  exp_matrix_binary <- exp_matrix * 1
  
  # Remove columns with species without data in a given experiment from the matrix
  keep <- sapply(colnames(exp_matrix_binary), function(x){
    exp_name %in% names(summaryAllSpAllExp[[x]])
  })
  exp_matrix_binary <- exp_matrix_binary[, which(keep)]
  
  # Remove rows not found in any species in a given experiment
  exp_matrix_binary <- exp_matrix_binary[which(rowSums(exp_matrix_binary) > 0), ]
  
  # Join the TLON and TMEL columns
  if(sum(colnames(exp_matrix_binary) %in% c("TLON", "TMEL")) > 0){
    if(sum(colnames(exp_matrix_binary) %in% c("TLON", "TMEL")) > 1){
      newVec <- rowSums(exp_matrix_binary[, which(colnames(exp_matrix_binary) %in% c("TLON", "TMEL"))]) > 0
      exp_matrix_binary <- cbind(exp_matrix_binary, newVec)
      colnames(exp_matrix_binary)[ncol(exp_matrix_binary)] <- "TETR"
      exp_matrix_binary <- exp_matrix_binary[, which(!colnames(exp_matrix_binary) %in% c("TLON", "TMEL"))]
    } else{
      colnames(exp_matrix_binary)[which(colnames(exp_matrix_binary) %in% c("TLON", "TMEL"))] <- "TETR"
    }
  }
  
  # Define PDF filename
  pdf_filename <- paste0("heatmap_unfiltered_", exp_name, ".pdf")
  
  # Create and save the heatmap
  pdf(pdf_filename, width = 8, height = 24)  # Adjust size as needed
  pheatmap(exp_matrix_binary,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           color = c("#ECEAE6", "#0a9396"),
           main = paste("Pathway Presence -", exp_name),
           legend=FALSE,
           #border_color = "grey",
           cellwidth = 8,
           cellheight = 8,
           fontsize = 8)
  dev.off()
}

```

Then, we want to check the 10 most represented pathways across species for each experiment 

[Per_Experiment_pathways_barplot.R](Pathways%20Functional%20Convergence%20203f4c09425a803f946dca91f0da54dc/Per_Experiment_pathways_barplot.r)

```r
# Set working directory
setwd("./Pathway_analysis")

# Load required libraries
library(dplyr)
library(stringr)
library(tibble)
library(pbmcapply)
library(Biostrings)
library(pheatmap)
library(KEGGREST)

################################################################################
# Step 1: Build KEGG pathway hierarchy
################################################################################

brite_raw <- keggGet("br:br08901")[[1]]
lines <- str_split(brite_raw, "\n")[[1]]

records <- list()
level1 <- level2 <- NULL

for (line in lines) {
  line <- str_trim(line)
  if (str_starts(line, "A")) {
    level1 <- str_remove(line, "^A\\s*")
  } else if (str_starts(line, "B")) {
    level2 <- str_remove(line, "^B\\s*")
  } else if (str_starts(line, "C")) {
    path_id <- str_extract(line, "\\d{5}")
    path_id <- paste0("path:map:", path_id)
    path_name <- str_trim(str_remove(line, "^C\\s*\\d{5}\\s*"))
    if (!is.na(path_id)) {
      records[[length(records) + 1]] <- tibble(
        Pathway_ID = path_id,
        Level3 = path_name,
        Level2 = level2,
        Level1 = level1
      )
    }
  }
}

# Combine and filter out non-relevant categories
pathway_df <- bind_rows(records)
pathway_df <- pathway_df[!pathway_df$Level1 %in% c("Human Diseases", "Drug Development"), ]

################################################################################
# Step 2: Load summarized annotation data
################################################################################

load("summaryAllSpAllExp.RData")

################################################################################
# Step 3: Build presence/absence matrices per experiment
################################################################################

# Pathway list
all_pathways <- sub("path:map:", "map", pathway_df$Pathway_ID)

# Keep all desired species (without collapsing TLON/TMEL)
desired_column_order <- c("CELE", "LMAR", "PEAU", "PLAE", "MISO", "ONUN", "SMED", "PTUR",
                          "SPEC", "TPIS", "PACU", "COCO", "EAND", "HMED", "LEPN", "TLON", "TMEL")

# List of experiments
exp_list <- c("CD", "CF", "Hyper", "Hypo", "OR", "UV0", "UV15D", "UV15L", "UV24D", "UV24L", "VL")

# Initialize list to hold all matrices
heatMaps_byExp_full <- list()

# Loop through experiments
for (exp_name in exp_list) {
  message("Processing: ", exp_name)
  
  countsbyExp <- lapply(summaryAllSpAllExp, function(df) {
    if (exp_name %in% colnames(df)) {
      relevant_rows <- df[(!df$pathways %in% c("", "-")) & df[[exp_name]], ]
      sapply(all_pathways, function(path_id) {
        sum(grepl(path_id, relevant_rows$pathways)) / sum(df[[exp_name]]) > 0.01
      })
    } else {
      NULL
    }
  })
  
  # Combine into matrix
  countsbyExp <- do.call(cbind, countsbyExp)
  colnames(countsbyExp) <- sub("SPEC1", "SPEC", colnames(countsbyExp))
  
  # Reorder columns to match desired order
  common_species <- intersect(desired_column_order, colnames(countsbyExp))
  countsbyExp <- countsbyExp[, common_species, drop = FALSE]
  rownames(countsbyExp) <- all_pathways
  if(sum(colnames(countsbyExp) %in% c("TLON", "TMEL")) > 0){
    if(sum(colnames(countsbyExp) %in% c("TLON", "TMEL")) > 1){
      newVec <- rowSums(countsbyExp[, which(colnames(countsbyExp) %in% c("TLON", "TMEL"))]) > 0
      countsbyExp <- cbind(countsbyExp, newVec)
      colnames(countsbyExp)[ncol(countsbyExp)] <- "TETR"
      countsbyExp <- countsbyExp[, which(!colnames(countsbyExp) %in% c("TLON", "TMEL"))]
    } else{
      colnames(countsbyExp)[which(colnames(countsbyExp) %in% c("TLON", "TMEL"))] <- "TETR"
    }
  }
  
  # Store in output list
  heatMaps_byExp_full[[exp_name]] <- countsbyExp
}

#Filter each experiment to keep pathways present in more than 50% of the species
filtered_heatMaps_byExp <- lapply(heatMaps_byExp_full, function(hm) {
  hm[rowMeans(hm) > 0.50, ]
})
filtered_heatMaps_byExp

#### Map the pathways IDs to the pathways names 
# Clean up the Pathway_ID in pathway_df to match the rownames in your heatmaps
pathway_df <- pathway_df %>%
  mutate(Pathway_ID_clean = gsub("path:map:", "map", Pathway_ID))

# Make a named vector to use for renaming
pathway_names <- setNames(pathway_df$Level3, pathway_df$Pathway_ID_clean)

# Apply renaming to each filtered heatmap
filtered_named_heatmaps <- lapply(filtered_heatMaps_byExp, function(hm) {
  rownames(hm) <- ifelse(rownames(hm) %in% names(pathway_names),
                         pathway_names[rownames(hm)],
                         rownames(hm))  # Keep original name if not found
  return(hm)
})

library(tidyverse)

# Create bar plots of top 10 most shared pathways per experiment
plots <- imap(filtered_named_heatmaps, function(mat, exp_name) {
  mat_numeric <- mat * 1  # Convert logicals to numeric
  
  tibble(Pathway = rownames(mat_numeric),
         Count = rowSums(mat_numeric)) %>%
    arrange(desc(Count)) %>%
    slice_head(n = 10) %>%
    ggplot(aes(x = reorder(Pathway, Count), y = Count)) +
    geom_col(fill = "#005f73") +
    coord_flip() +
    labs(title = paste("Top Shared Pathways in", exp_name),
         x = "Pathway",
         y = "Number of Species") +
    theme_minimal()
})

library(patchwork)

# Combine all plots side by side (adjust ncol to control layout)
combined_plot <- wrap_plots(plots, ncol = length(plots)/2)  # All in one row

# Display
print(combined_plot)

# Save as a wide figure
ggsave("top_pathways_side_by_side.png", combined_plot, width = 15, height = 5)

library(tidyverse)
library(ggpubr)

# Combine all top 10 pathways per experiment into a single data frame
all_data <- imap_dfr(filtered_named_heatmaps, function(mat, exp_name) {
  tibble(
    Pathway = rownames(mat),
    Count = rowSums(mat * 1),
    Experiment = exp_name
  )
})

# Select top 10 pathways per experiment
top_pathways <- all_data %>%
  group_by(Experiment) %>%
  slice_max(order_by = Count, n = 10, with_ties = FALSE) %>%
  ungroup()
top_pathways$Pathway <- factor(top_pathways$Pathway, levels = rev(unique(top_pathways$Pathway)))

# Plot them together in one figure with facets
allPlots <- lapply(unique(top_pathways$Experiment), function(x){
  temp_pathways <- top_pathways[top_pathways$Experiment == x, ]
  temp_pathways$Pathway <- factor(temp_pathways$Pathway, levels = rev(unique(temp_pathways$Pathway)))
  out <- ggplot(temp_pathways, aes(x = Pathway, y = Count)) +
    geom_col(fill = "#005f73", width = 0.6) +
    coord_flip() +
    scale_y_continuous(breaks = seq(0, 16, by = 2), limits = c(0, 16)) +
    labs(
      title = x,
      x = "Pathway",
      y = "Number of Species"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )
})
final_plot <- ggarrange(plotlist = allPlots, common.legend = TRUE, nrow = 3, ncol = 4)
annotate_figure(final_plot, top = text_grob("Top 10 Shared Pathways per Experiment", face = "bold", size = 14))

```