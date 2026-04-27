# Pathway Functional Convergence (All Genes Analysis)

To assess functional convergence at the pathway level, KEGG pathway annotations were assigned to all genes belonging to significant WGCNA modules across species and stress experiments. KEGG annotations were obtained using eggNOG-mapper v2.1.12 and integrated with gene–protein mappings derived from the co-expression network analysis.

For each species and experiment, pathway representation was quantified using all significant module-associated genes. Pathways were considered relevant if they accounted for at least 1% of the genes within a given species–experiment context.

Pathway representation is expressed as percentages, allowing direct comparison of pathway contribution across species and experimental conditions.

We have summarised the results of eggNOG-mapper in this file: ids_pathways_with_exp.RData

First, we extracted the pathways that are present in more than 80% of the species, combining all experiments and visualised the percentage of present heatmap each time, plus per experiment using the same set of pathways (Figure 4a, Fig S1). Then we calculated the heatmap percentage  present per experiment without filtering.

First, for each species and experiment, we calculated the percentage of significant module-associated genes assigned to each KEGG pathway. A pathway was considered present in a species–experiment context when it represented more than 1% of the genes. We then combined all experiments and retained pathways present in more than 80% of species. The corresponding pathway percentages were visualized in a combined heatmap, and the same filtered pathway set was used to generate per-experiment heatmaps (Fig. 4a, Fig. S129a).

In addition, we generated unfiltered per-experiment heatmaps, showing the percentage representation of all KEGG pathways detected in each experiment (Figs. S130-S140)

```r
# Set working directory
setwd("/Pathway_analysis_ALL_GENES")

# Load required libraries
library(dplyr)
library(stringr)
library(tibble)
library(pbmcapply)
library(Biostrings)
library(pheatmap)
library(KEGGREST)

output_dir <- file.path(getwd(), "Pathway_Percentage_Heatmaps_white")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Build KEGG pathway hierarchy
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

# Load summarized annotation data
load("ids_pathways_with_exp.RData")

# Step 3: Build percentage matrices per experiment
# Pathway list
all_pathways <- sub("path:map:", "map", pathway_df$Pathway_ID)

# Keep all desired species
desired_column_order <- c("CELE", "LMAR", "PEAU", "PLAE", "MISO", "ONUN", "SMED", "PTUR",
                          "SPEC", "TPIS", "PACU", "COCO", "EAND", "HMED", "LEPN", "TLON", "TMEL")

# List of experiments
exp_list <- c("CD", "CF", "Hyper", "Hypo", "OR", "UV0", "UV15D", "UV15L", "UV24D", "UV24L", "VL")

# Initialize lists
heatMaps_byExp_percent <- list()

for (exp_name in exp_list) {
  message("Processing: ", exp_name)
  
  percent_byExp <- lapply(ids_pathways_with_exp, function(df) {
    if (exp_name %in% colnames(df)) {
      relevant_rows <- df[(!df$pathways %in% c("", "-")) & df[[exp_name]], ]
      total_genes <- sum(df[[exp_name]], na.rm = TRUE)
      
      sapply(all_pathways, function(path_id) {
        if (total_genes == 0) {
          NA_real_
        } else {
          100 * sum(grepl(path_id, relevant_rows$pathways)) / total_genes
        }
      })
    } else {
      rep(NA_real_, length(all_pathways))
    }
  })
  
  percent_byExp <- do.call(cbind, percent_byExp)
  colnames(percent_byExp) <- sub("SPEC1", "SPEC", names(ids_pathways_with_exp))
  common_species <- intersect(desired_column_order, colnames(percent_byExp))
  percent_byExp <- percent_byExp[, common_species, drop = FALSE]
  rownames(percent_byExp) <- all_pathways
 
  heatMaps_byExp_percent[[exp_name]] <- percent_byExp
}

# Save percentage matrices

for (exp in names(heatMaps_byExp_percent)) {
  write.csv(
    heatMaps_byExp_percent[[exp]],
    file.path(output_dir, paste0("PathwayPercent_", exp, ".csv")),
    row.names = TRUE
  )
}

# Combine all experiment matrices into one big matrix

# Percentage combined matrix
combinedPercentMatrix <- do.call(cbind, lapply(names(heatMaps_byExp_percent), function(exp_name) {
  mat <- heatMaps_byExp_percent[[exp_name]]
  colnames(mat) <- paste(colnames(mat), exp_name, sep = "_")
  return(mat)
}))

# Merge experiments by species
# For percentages: merge by OR translated as MAX percentage across experiments
# TLON + TMEL are merged into TETR

all_colnames <- colnames(combinedPercentMatrix)
base_species <- unique(sub("_.*", "", all_colnames))
base_species <- unique(sub("TLON|TMEL", "TETR", base_species))

# Merged percentage matrix
mergedPercentMatrix <- sapply(base_species, function(sp) {
  if (sp == "TETR") {
    cols <- grep("^(TLON|TMEL)_", all_colnames)
  } else {
    cols <- grep(paste0("^", sp, "_"), all_colnames)
  }
  
  if (length(cols) == 0) {
    rep(NA_real_, nrow(combinedPercentMatrix))
  } else {
    apply(combinedPercentMatrix[, cols, drop = FALSE], 1, max, na.rm = TRUE)
  }
})

mergedPercentMatrix <- as.matrix(mergedPercentMatrix)
rownames(mergedPercentMatrix) <- rownames(combinedPercentMatrix)

# Merged percent matrix
write.csv(
  mergedPercentMatrix,
  file.path(output_dir, "PathwayPercent_MergedBySpecies_TETR.csv"),
  row.names = TRUE
)

# Filter pathways by keeping pathways present in >80% of species

threshold <- 0.8
species_count <- ncol(mergedPercentMatrix)
pathway_freq <- rowSums(mergedPercentMatrix, na.rm = TRUE) / species_count

filteredPercentMatrix <- mergedPercentMatrix[pathway_freq > threshold, , drop = FALSE]

write.csv(
  filteredPercentMatrix,
  file.path(output_dir, "PathwayPercent_Filtered_1%_80.csv"),
  row.names = TRUE
)

# Keep only filtered pathways in each experiment
filtered_pathway_ids <- rownames(filteredPercentMatrix)

heatMaps_byExp_filtered_percent <- list()

for (exp_name in names(heatMaps_byExp_percent)) {
  exp_matrix <- heatMaps_byExp_percent[[exp_name]]
  pathway_ids_in_exp <- rownames(exp_matrix)
  matching_pathway_ids <- intersect(filtered_pathway_ids, pathway_ids_in_exp)
  exp_matrix_filtered <- exp_matrix[matching_pathway_ids, , drop = FALSE]
  heatMaps_byExp_filtered_percent[[exp_name]] <- exp_matrix_filtered
}

# Function to map pathway IDs to full pathway names
map_full_names <- function(pathway_ids) {
  full_names <- sapply(pathway_ids, function(pid) {
    match_idx <- which(sub("path:map:", "map", pathway_df$Pathway_ID) == pid)
    if (length(match_idx) > 0) {
      pathway_df$Level3[match_idx]
    } else {
      pid
    }
  })
  return(full_names)
}

# Plot merged filtered heatmap with full pathway names
# AND actual number of pathways per species in column labels

pathway_ids <- rownames(filteredPercentMatrix)
full_names <- map_full_names(pathway_ids)

filteredMatrix_named <- filteredPercentMatrix

# Remove only the pathway called "Metabolic pathways"
metabolic_ids <- pathway_df$Pathway_ID[pathway_df$Level3 == "Metabolic pathways"]
metabolic_ids <- sub("path:map:", "map", metabolic_ids)

filteredMatrix_named <- filteredMatrix_named[
  !rownames(filteredMatrix_named) %in% metabolic_ids,
  , drop = FALSE
]

# NOW get pathway IDs after filtering
pathway_ids <- rownames(filteredMatrix_named)
full_names <- map_full_names(pathway_ids)

rownames(filteredMatrix_named) <- full_names

# Check dimensions before plotting
print(dim(filteredMatrix_named))
print(colnames(filteredMatrix_named))
print(head(filteredMatrix_named[, , drop = FALSE]))

# Stop early if matrix is empty
if (nrow(filteredMatrix_named) == 0 || ncol(filteredMatrix_named) == 0) {
  stop("filteredMatrix_named is empty: no pathways passed the filter or no species columns are present.")
}

# Count actual number of filtered pathways present per species
# Presence follows the 1% rule
species_pathway_counts <- colSums(filteredPercentMatrix > 1, na.rm = TRUE)
print(species_pathway_counts)

# Make sure counts align with columns
species_pathway_counts <- species_pathway_counts[colnames(filteredMatrix_named)]

# Add counts directly to column labels
new_colnames <- paste0(colnames(filteredMatrix_named), " (n=", species_pathway_counts, ")")
print(new_colnames)
colnames(filteredMatrix_named) <- new_colnames

# Use the requested common visualization scale

max_pct <- ceiling(max(filteredMatrix_named, na.rm = TRUE))

pdf(file.path(output_dir, "heatmap_allExperiments_percent.pdf"), width = 10, height = 10)

pheatmap(
  filteredMatrix_named,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white",
    "#EDC8C5",
    "#F0A562",
    "#A54055",
    "#6B2347",
    "#4C1B6B",
    "#1A1E50"
  ))(100),
  breaks = seq(0, max_pct, length.out = 101),
  cellwidth = 10,
  cellheight = 10,
  fontsize = 6,
  main = "Pathways in >80% of Species (>1%)",
  silent = FALSE, angle_col = 90
)

dev.off()

# Plot filtered heatmaps for each experiment
# Add actual number of pathways per species to column labels

first_plot <- TRUE

for (exp_name in names(heatMaps_byExp_filtered_percent)) {
  exp_matrix_filtered <- heatMaps_byExp_filtered_percent[[exp_name]]
  
  # Remove only "Metabolic pathways"
  metabolic_ids <- pathway_df$Pathway_ID[pathway_df$Level3 == "Metabolic pathways"]
  metabolic_ids <- sub("path:map:", "map", metabolic_ids)
  
  exp_matrix_filtered <- exp_matrix_filtered[
    !rownames(exp_matrix_filtered) %in% metabolic_ids,
    , drop = FALSE
  ]
  
  # Map names AFTER filtering
  pathway_ids <- rownames(exp_matrix_filtered)
  full_names <- map_full_names(pathway_ids)
  rownames(exp_matrix_filtered) <- full_names
  
  pathway_ids <- rownames(exp_matrix_filtered)
  full_names <- map_full_names(pathway_ids)
  rownames(exp_matrix_filtered) <- full_names
  
  # Remove columns with species without data in a given experiment
  keep <- sapply(colnames(exp_matrix_filtered), function(x) {
    exp_name %in% names(ids_pathways_with_exp[[x]])
  })
  exp_matrix_filtered <- exp_matrix_filtered[, which(keep), drop = FALSE]
  
  # Join TLON and TMEL columns as TETR using max percentage
  if (sum(colnames(exp_matrix_filtered) %in% c("TLON", "TMEL")) > 0) {
    if (sum(colnames(exp_matrix_filtered) %in% c("TLON", "TMEL")) > 1) {
      tetr_vec <- apply(
        exp_matrix_filtered[, colnames(exp_matrix_filtered) %in% c("TLON", "TMEL"), drop = FALSE],
        1,
        function(x) {
          if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
        }
      )
      exp_matrix_filtered <- cbind(exp_matrix_filtered, TETR = tetr_vec)
      exp_matrix_filtered <- exp_matrix_filtered[, !colnames(exp_matrix_filtered) %in% c("TLON", "TMEL"), drop = FALSE]
    } else {
      colnames(exp_matrix_filtered)[colnames(exp_matrix_filtered) %in% c("TLON", "TMEL")] <- "TETR"
    }
  }
  
  # Clean matrix
  exp_matrix_filtered <- as.matrix(exp_matrix_filtered)
  mode(exp_matrix_filtered) <- "numeric"
  exp_matrix_filtered[!is.finite(exp_matrix_filtered)] <- NA_real_
  rownames(exp_matrix_filtered) <- make.unique(rownames(exp_matrix_filtered))
  
  # Count actual number of pathways present per species (>1%)
  species_pathway_counts <- colSums(exp_matrix_filtered > 1, na.rm = TRUE)
  
  # Put counts on the same line, shorter label
  colnames(exp_matrix_filtered) <- paste0(
    colnames(exp_matrix_filtered),
    " (", species_pathway_counts[colnames(exp_matrix_filtered)], ")"
  )
  
  max_pct <- ceiling(max(exp_matrix_filtered, na.rm = TRUE))
  pdf_filename <- paste0("heatmap_filtered_", exp_name, "_percent.pdf")
  pdf(file.path(output_dir, pdf_filename), width = 8, height = 8)
  
  pheatmap(
    exp_matrix_filtered,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = first_plot,
    show_colnames = TRUE,
    color = colorRampPalette(c("white",
      "#EDC8C5",
      "#F0A562",
      "#A54055",
      "#6B2347",
      "#4C1B6B",
      "#1A1E50"
    ))(100),
    breaks = seq(0, max_pct, length.out = 101),
    main = paste("Pathway Percentage -", exp_name),
    legend = TRUE,
    cellwidth = 9,
    cellheight = 8,
    fontsize = 8,
    fontsize_col = 6,
    angle_col = 90
  )
  
  dev.off()
  first_plot <- FALSE
}

################################################################################
# Plot unfiltered heatmaps for each experiment
# Add actual number of pathways per species to column labels
################################################################################

for (exp_name in names(heatMaps_byExp_percent)) {
  exp_matrix <- heatMaps_byExp_percent[[exp_name]]
  
  pathway_ids <- rownames(exp_matrix)
  full_names <- map_full_names(pathway_ids)
  rownames(exp_matrix) <- full_names
  
  # Remove columns with species without data in a given experiment
  keep <- sapply(colnames(exp_matrix), function(x) {
    exp_name %in% names(ids_pathways_with_exp[[x]])
  })
  exp_matrix <- exp_matrix[, which(keep), drop = FALSE]
  
  # Remove rows not found in any species in a given experiment
  exp_matrix <- exp_matrix[which(rowSums(exp_matrix, na.rm = TRUE) > 0), , drop = FALSE]
  
  # Join TLON and TMEL columns as TETR using max percentage
  if (sum(colnames(exp_matrix) %in% c("TLON", "TMEL")) > 0) {
    if (sum(colnames(exp_matrix) %in% c("TLON", "TMEL")) > 1) {
      tetr_vec <- apply(
        exp_matrix[, colnames(exp_matrix) %in% c("TLON", "TMEL"), drop = FALSE],
        1,
        function(x) {
          if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
        }
      )
      exp_matrix <- cbind(exp_matrix, TETR = tetr_vec)
      exp_matrix <- exp_matrix[, !colnames(exp_matrix) %in% c("TLON", "TMEL"), drop = FALSE]
    } else {
      colnames(exp_matrix)[colnames(exp_matrix) %in% c("TLON", "TMEL")] <- "TETR"
    }
  }
  
  # Clean matrix
  exp_matrix <- as.matrix(exp_matrix)
  mode(exp_matrix) <- "numeric"
  exp_matrix[!is.finite(exp_matrix)] <- NA_real_
  rownames(exp_matrix) <- make.unique(rownames(exp_matrix))
  
  # Count actual number of pathways present per species (>1%)
  species_pathway_counts <- colSums(exp_matrix > 1, na.rm = TRUE)
  
  # Add counts on the SAME line (like filtered)
  colnames(exp_matrix) <- paste0(
    colnames(exp_matrix),
    " (", species_pathway_counts[colnames(exp_matrix)], ")"
  )
  
  
  n_rows <- nrow(exp_matrix)
  max_pct <- ceiling(max(exp_matrix, na.rm = TRUE))
  
  # Set height proportional to number of pathways
  pdf_height <- max(24, n_rows * 0.12)
  pdf_filename <- paste0("heatmap_unfiltered_", exp_name, "_percent.pdf")
  pdf(file.path(output_dir, pdf_filename), width = 8, height = pdf_height)  
  pheatmap(
    exp_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    color = colorRampPalette(c("white",
      "#EDC8C5",
      "#F0A562",
      "#A54055",
      "#6B2347",
      "#4C1B6B",
      "#1A1E50"
    ))(100),
    breaks = seq(0, max_pct, length.out = 101),
    main = paste("Pathway Percentage -", exp_name),
    legend = TRUE,
    cellwidth = 8,
    cellheight = 8,
    fontsize = 6,
    angle_col = 90
  )
  
  dev.off()
}

################################################################################
# Done
################################################################################

message("All percentage matrices and heatmaps were saved in: ", output_dir)
```

Then, we want to check the 10 most represented pathways across species for each experiment (Fig. 4b, Fig. S129b).

```r
################################################################################
# Setup
################################################################################
setwd("Pathway_analysis_ALL_GENES/Pathway_Percentage_Heatmaps")

library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(purrr)
library(KEGGREST)
library(ggplot2)
library(ggpubr)

# Build KEGG pathway hierarchy (BRITE)
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

pathway_df <- bind_rows(records) %>%
  filter(!Level1 %in% c("Human Diseases", "Drug Development")) %>%
  mutate(Pathway_ID_clean = gsub("path:map:", "map", Pathway_ID))

# Named vector: mapXXXXX -> pathway name
pathway_names <- setNames(pathway_df$Level3, pathway_df$Pathway_ID_clean)

# Load summarized annotation data
load("./ids_pathways_with_exp.RData")

# Pathway IDs used for rownames in matrices
all_pathways <- pathway_df$Pathway_ID_clean

# Desired species order 
desired_column_order <- c(
  "CELE", "LMAR", "PEAU", "PLAE", "MISO", "ONUN", "SMED", "PTUR",
  "SPEC", "TPIS", "PACU", "COCO", "EAND", "HMED", "LEPN", "TLON", "TMEL"
)

# Experiments
exp_list <- c("CD", "CF", "Hyper", "Hypo", "OR", "UV0", "UV15D", "UV15L", "UV24D", "UV24L", "VL")

# Collapse TLON/TMEL into TETR
collapse_tetr <- function(mat) {
  if (!any(colnames(mat) %in% c("TLON", "TMEL"))) return(mat)
  
  if (sum(colnames(mat) %in% c("TLON", "TMEL")) > 1) {
    tetr <- rowSums(mat[, colnames(mat) %in% c("TLON", "TMEL"), drop = FALSE]) > 0
    mat <- cbind(mat, TETR = tetr)
    mat <- mat[, !colnames(mat) %in% c("TLON", "TMEL"), drop = FALSE]
  } else {
    colnames(mat)[colnames(mat) %in% c("TLON", "TMEL")] <- "TETR"
  }
  
  mat
}

heatmaps_by_exp <- list()

for (exp_name in exp_list) {
  message("Processing: ", exp_name)
  
  per_species <- lapply(ids_pathways_with_exp, function(df) {
    if (!exp_name %in% colnames(df)) return(NULL)
    
    # rows that belong to this experiment and have pathway annotation
    relevant_rows <- df[(!df$pathways %in% c("", "-")) & df[[exp_name]], , drop = FALSE]
    
    denom <- sum(df[[exp_name]], na.rm = TRUE)
    if (is.na(denom) || denom == 0) {
      # if no genes flagged for that experiment in this species, return FALSE for all pathways
      return(setNames(rep(FALSE, length(all_pathways)), all_pathways))
    }
    
    # present if >1% of genes in that species/experiment have that pathway
    out <- sapply(all_pathways, function(path_id) {
      sum(grepl(path_id, relevant_rows$pathways)) / denom > 0.01
    })
    
    out
  })
  
  mat <- do.call(cbind, per_species)
  colnames(mat) <- sub("SPEC1", "SPEC", colnames(mat))
  
  # reorder columns
  common_species <- intersect(desired_column_order, colnames(mat))
  mat <- mat[, common_species, drop = FALSE]
  
  rownames(mat) <- all_pathways
  
  # collapse TLON/TMEL -> TETR
  mat <- collapse_tetr(mat)
  
  heatmaps_by_exp[[exp_name]] <- mat
}

# Filter pathways present in >50% species (per experiment)
filtered_heatmaps <- lapply(heatmaps_by_exp, function(hm) {
  hm[rowMeans(hm) > 0.50, , drop = FALSE]
})

# Rename pathway IDs to pathway names (Level3)
filtered_named_heatmaps <- lapply(filtered_heatmaps, function(hm) {
  rn <- rownames(hm)
  rownames(hm) <- ifelse(rn %in% names(pathway_names), pathway_names[rn], rn)
  hm
})

# Select top 10 per experiment
all_data <- imap_dfr(filtered_named_heatmaps, function(mat, exp_name) {
  tibble(
    Experiment = exp_name,
    Pathway    = rownames(mat),
    Count      = rowSums(mat * 1)
  )
})

top_pathways <- all_data %>%
  group_by(Experiment) %>%
  slice_max(order_by = Count, n = 10, with_ties = FALSE) %>%
  ungroup()

# Plot (multi-panel, ordered within each experiment)
make_one_plot <- function(df_exp) {
  df_exp <- df_exp %>%
    arrange(Count) %>%
    mutate(Pathway = factor(Pathway, levels = Pathway))  # order within panel
  
  ggplot(df_exp, aes(x = Pathway, y = Count)) +
    geom_col(fill = "#C96A74", width = 0.6) +
    coord_flip() +
    scale_y_continuous(breaks = seq(0, 16, by = 2), limits = c(0, 16)) +
    labs(
      title = unique(df_exp$Experiment),
      x = "Pathway",
      y = "Number of Species"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 8)
    )
}

plot_list <- split(top_pathways, top_pathways$Experiment) %>%
  lapply(make_one_plot)

final_plot <- ggarrange(plotlist = plot_list, nrow = 3, ncol = 4)
final_plot <- annotate_figure(
  final_plot,
  top = text_grob("Top 10 Shared Pathways per Experiment", face = "bold", size = 14)
)

print(final_plot)

# Save
ggsave("top_10_shared_pathways_per_experiment.png", final_plot, width = 16, height = 10, dpi = 300)

ggsave(
  "top_10_shared_pathways_per_experiment.pdf",
  final_plot,
  width = 25,
  height = 10,
)

```

### Module composition statistics based on KEGG pathway annotations

To describe the functional composition of significant WGCNA modules, we calculated KEGG pathway annotation metrics for each module across all species. For each significant module, gene identifiers were matched to translated protein identifiers and their corresponding KEGG pathway annotations.

For every species-specific module, we recorded the total number of transcripts in the module, the number of translated transcripts, the number of proteins with KEGG annotations, the number of unique KEGG pathways represented, and pathway richness. Pathway richness was calculated as the number of unique KEGG pathways divided by the total number of transcripts in the module. (Table S58).

```python
from pathlib import Path
import pandas as pd
import re

base_dir = Path("/path/to/base")

work_dir = Path("/Pathway_analysis_ALL_GENES")

out_dir = work_dir / "Module_KEGG_Composition"
out_dir.mkdir(parents=True, exist_ok=True)

# READ KEGG ANNOTATION

prot2paths = work_dir / "ids_pathways_with_exp_flat.csv"

prot2paths2genes = pd.read_csv(prot2paths, usecols=[0, 1])
prot2paths2genes.columns = ["protein_id", "pathway"]

# Derive gene ID from protein ID
prot2paths2genes["gene_id"] = prot2paths2genes["protein_id"].str.replace(
    r"\.p\d+$", "", regex=True
)

# SPECIES METADATA
species_meta = [
    ("Annelida",        "Arenicola_marina",            "COCO"),
    ("Annelida",        "Eisenia_andrei",              "EAND"),
    ("Annelida",        "Hirudo_medicinalis",          "HMED"),
    ("Arthropoda",      "Ligia_oceanica",              "MISO"),
    ("Arthropoda",      "Porcellio_laevis",            "PLAE"),
    ("Mollusca",        "Phorcus_turbinatus",          "PTUR"),
    ("Mollusca",        "Physella_acuta",              "PACU"),
    ("Mollusca",        "Siphonaria_pectinata",        "SPEC"),
    ("Mollusca",        "Theba_pisana",                "TPIS"),
    ("Nematoda",        "Litoditis_marina",            "LMAR"),
    ("Nematoda",        "Caenorhabditis_elegans",      "CELE"),
    ("Nemertea",        "Tetrastemma_melanocephalum",  "TMEL"),
    ("Nemertea",        "Tetrastemma_longissimum",     "TLON"),
    ("Nemertea",        "Leptonemertes_chalicophora",  "LEPN"),
    ("Onychophora",     "Peripatoides_aurorbis",       "PEAU"),
    ("Platyhelminthes", "Obama_nungara",               "ONUN"),
    ("Platyhelminthes", "Schmidtea_mediterranea",      "SMED"),
]

# FIND MODULE FILES
pattern = re.compile(
    r"^(?P<code>[^_]+)_(?P<module>.+?)_significant_module_geneList_for_(?P<condition>.+)\.txt$"
)

records = []

for phylum, species_folder, code in species_meta:
    mod_dir = (
        base_dir
        / phylum
        / species_folder
        / "WGCNA"
        / "3.Relate_modules_to_external_traits"
        / "1.signifMods2genes"
    )

    if not mod_dir.exists():
        print(f"WARNING: directory not found: {mod_dir}")
        continue

    for fp in mod_dir.glob("*.txt"):
        m = pattern.match(fp.name)
        if not m:
            continue

        records.append({
            "phylum": phylum,
            "species_folder": species_folder,
            "species_code": code,
            "module_color": m.group("module"),
            "condition": m.group("condition"),
            "filepath": str(fp),
        })

files_df = pd.DataFrame(records)

# Keep one representative file per species/module.
# This avoids counting the same module more than once if it is significant
# in multiple experimental conditions.
representative_files_df = (
    files_df
    .sort_values("condition")
    .drop_duplicates(subset=["species_code", "module_color"])
    .reset_index(drop=True)
)

print("Total unique modules:", len(representative_files_df))

def expand_module_pathways(module_df):
    """Split comma-separated KEGG pathway annotations into one row per pathway."""
    df = module_df.copy()

    df["pathway"] = df["pathway"].replace("-", pd.NA)
    df["pathway"] = df["pathway"].astype("string")
    df["pathway"] = df["pathway"].str.split(",")

    df = df.explode("pathway", ignore_index=True)
    df["pathway"] = df["pathway"].str.strip()

    return df

def summarize_one_module(module_df):
    """Calculate KEGG composition statistics for one module."""
    expanded = expand_module_pathways(module_df)

    n_genes_total = module_df["gene_id"].nunique()

    n_genes_translated = module_df.loc[
        module_df["protein_id"].notna(),
        "gene_id"
    ].nunique()

    n_annotated_proteins = expanded.loc[
        expanded["pathway"].notna(),
        "protein_id"
    ].nunique()

    n_unique_pathways = expanded.loc[
        expanded["pathway"].notna(),
        "pathway"
    ].nunique()

    pathway_richness = (
        n_unique_pathways / n_genes_total
        if n_genes_total > 0
        else pd.NA
    )

    return {
        "Phylum": module_df["phylum"].iloc[0],
        "Species": module_df["species_folder"].iloc[0],
        "Species_code": module_df["species_code"].iloc[0],
        "Module color": module_df["module_color"].iloc[0],
        "Number of transcripts in module": n_genes_total,
        "Number of proteins (translated transcripts)": n_genes_translated,
        "Annotated proteins": n_annotated_proteins,
        "Number of unique pathways": n_unique_pathways,
        "Pathway richness": pathway_richness,
    }

# SUMMARIZE MODULES
all_module_summaries = []

for _, row in representative_files_df.iterrows():
    fp = Path(row["filepath"])

    try:
        with open(fp, "r") as f:
            genes = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"WARNING: could not read {fp}: {e}")
        continue

    module_df = pd.DataFrame({"gene_id": genes})

    # Merge gene IDs with protein IDs and KEGG pathways
    module_df = module_df.merge(
        prot2paths2genes[["gene_id", "protein_id", "pathway"]],
        on="gene_id",
        how="left"
    )

    module_df.insert(0, "module_color", row["module_color"])
    module_df.insert(0, "species_code", row["species_code"])
    module_df.insert(0, "species_folder", row["species_folder"])
    module_df.insert(0, "phylum", row["phylum"])

    all_module_summaries.append(summarize_one_module(module_df))

all_module_summary_df = pd.DataFrame(all_module_summaries)

output_file = out_dir / "all_module_composition_summary.csv"
all_module_summary_df.to_csv(output_file, index=False)
```