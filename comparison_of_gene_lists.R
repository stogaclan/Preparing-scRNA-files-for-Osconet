library(VennDiagram)

# Pools the gene names from all per-cluster *Comm.csv files produced by
# 4-CommunityExtract.R, separately for two conditions, then compares them.
#
# Outputs three gene lists (A only, B only, shared), a summary table, and
# a Venn diagram.


####################################################################
##                          EDIT SECTION                          ##
####################################################################

# folder holding the *Comm.csv files
resultsDir <- "~/OscoNet-Master-copy/Cavone2021_completed_ollies_ubuntu/Results"

# filename prefix for each condition, up to and including the 'D' that
# precedes the cluster number. Including the 'D' matters: a prefix of
# "Cavone2021_her4_les" alone would be fine here, but for cases like
# 'WT' and 'WTles' the shorter name is a substring of the longer one.
prefixA <- "Cavone2021_her4_naiveD"
prefixB <- "Cavone2021_her4_lesD"

# labels used in the Venn and in the output filenames
labelA <- "naive"
labelB <- "lesioned"

# where to write the results
outDir <- "~/OscoNet-Master-copy/Cavone2021_completed_ollies_ubuntu/Comparison"

# column holding the gene names
geneCol <- "GeneName"

####################################################################
##                  Do not modify after this line                 ##
####################################################################

resultsDir <- path.expand(resultsDir)
outDir     <- path.expand(outDir)

if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

# ---- collect genes per condition ---------------------------------
# returns a data.frame of gene / cluster so we keep track of which
# cluster each gene came from, not just the pooled set

collect <- function(prefix, label) {

  files <- list.files(resultsDir,
                      pattern = paste0("^", prefix, "[0-9]+_.*Comm\\.csv$"),
                      full.names = TRUE)

  if (length(files) == 0) {
    stop("No files matched prefix '", prefix, "' in ", resultsDir)
  }

  message("\n", label, ": ", length(files), " cluster files")

  out <- do.call(rbind, lapply(files, function(f) {
    d <- read.csv(f, check.names = FALSE)
    if (!geneCol %in% colnames(d)) {
      stop("Column '", geneCol, "' not found in ", basename(f),
           " - columns are: ", paste(colnames(d), collapse = ", "))
    }
    cluster <- sub(paste0("^", prefix), "", basename(f))
    cluster <- sub("_.*$", "", cluster)
    message("  ", basename(f), "  genes: ", nrow(d))
    data.frame(gene = as.character(d[[geneCol]]),
               cluster = cluster,
               stringsAsFactors = FALSE)
  }))

  message("  total rows: ", nrow(out),
          "   unique genes: ", length(unique(out$gene)))
  out
}

allA <- collect(prefixA, labelA)
allB <- collect(prefixB, labelB)

genesA <- sort(unique(allA$gene))
genesB <- sort(unique(allB$gene))

# ---- compare -----------------------------------------------------

onlyA  <- setdiff(genesA, genesB)
onlyB  <- setdiff(genesB, genesA)
shared <- intersect(genesA, genesB)

message("\n--- comparison ---")
message(labelA, " total:  ", length(genesA))
message(labelB, " total:  ", length(genesB))
message(labelA, " only:   ", length(onlyA))
message(labelB, " only:   ", length(onlyB))
message("shared:        ", length(shared))

# ---- write gene lists --------------------------------------------

writeLines(onlyA,  file.path(outDir, paste0(labelA, "_only.txt")))
writeLines(onlyB,  file.path(outDir, paste0(labelB, "_only.txt")))
writeLines(shared, file.path(outDir, "shared.txt"))

# ---- summary table -----------------------------------------------
# for each gene: which set it falls in, and how many clusters of each
# condition it appeared in

nClusA <- table(unique(allA[, c("gene", "cluster")])$gene)
nClusB <- table(unique(allB[, c("gene", "cluster")])$gene)

allGenes <- sort(union(genesA, genesB))

summaryTab <- data.frame(
  gene = allGenes,
  set = ifelse(allGenes %in% shared, "shared",
        ifelse(allGenes %in% onlyA, paste0(labelA, "_only"),
                                    paste0(labelB, "_only"))),
  clusters_A = as.integer(nClusA[allGenes]),
  clusters_B = as.integer(nClusB[allGenes]),
  stringsAsFactors = FALSE
)
summaryTab[is.na(summaryTab)] <- 0
colnames(summaryTab)[3:4] <- paste0("clusters_", c(labelA, labelB))

write.csv(summaryTab, file.path(outDir, "gene_comparison.csv"),
          row.names = FALSE)

# ---- Venn diagram ------------------------------------------------

vennFile <- file.path(outDir, paste0("venn_", labelA, "_vs_", labelB, ".png"))

futile.logger::flog.threshold(futile.logger::ERROR,
                              name = "VennDiagramLogger")

venn.diagram(
  x = setNames(list(genesA, genesB), c(labelA, labelB)),
  filename = vennFile,
  imagetype = "png",
  height = 1800, width = 2200, resolution = 300,
  # fixed, equal circles: by default venn.diagram scales each circle to its
  # set size and shifts the overlap to be area-proportional, which makes
  # diagrams from different datasets hard to compare side by side
  scaled = FALSE,
  euler.d = FALSE,
  fill = c("#7fb3d5", "#f1948a"),
  alpha = 0.6,
  lwd = 1,
  cex = 1.4,
  cat.cex = 1.4,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  margin = 0.1
)

message("\nWritten to ", outDir)
message("  ", labelA, "_only.txt   (", length(onlyA), " genes)")
message("  ", labelB, "_only.txt   (", length(onlyB), " genes)")
message("  shared.txt        (", length(shared), " genes)")
message("  gene_comparison.csv")
message("  ", basename(vennFile))
