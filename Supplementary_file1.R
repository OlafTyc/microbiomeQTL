#!/usr/bin/env Rscript
options(error=function()traceback(2))
#args = commandArgs(trailingOnly=TRUE)
#inp_folder <- args[[1]]
#prefix <- args[[2]]
inp_folder <- "C:/Users/Lokhorst/Downloads"
prefix <- "tomato3"

start.time <- Sys.time()

library(devtools)
#install.packages("qtl", repos = "http://cran.us.r-project.org")
#install.packages("qtl2", repos = "https://rqtl.org/qtl2cran")
library(qtl)
library(qtl2)
library(parallel)
library(MASS)
library(base)
#install.packages('RColorBrewer')
library(RColorBrewer)
library(plyr)

# Load OTU table
otu_table <- read.table(file.path(inp_folder, "counts_more_samples.txt"), sep = "\t", comment.char="", header=T, check.names=F)
colnames(otu_table) <- tolower(colnames(otu_table))
colnames(otu_table)[[1]] <- "id"

# Get the taxonomy
tax <- otu_table[,c(1, NCOL(otu_table))]

# Remove OTUs that are not likely to be of microbial origin
toMatch <- c("k__Unclassified", "f__Mitochondria", "o__Chloroplast")
remove_otus <- grep(paste(toMatch, collapse="|"), otu_table[,NCOL(otu_table)], fixed=F, value=F)
otu_table_removed <- otu_table[-remove_otus,]
tax1 <- tax[-remove_otus,]

# First pull out the abundances alone and convert to as.numeric (they were stored in as factors in a dataframe)
clean_otu <- otu_table_removed[1:NROW(otu_table_removed),2:(NCOL(otu_table_removed)-1)]
clean_otus <- data.matrix(clean_otu, rownames.force = NA)
# calculate sums
sums <- apply(clean_otus, 2, sum)
pic_file <- file.path(inp_folder, sprintf("%s_OTU_distribution.jpeg", prefix))
jpeg(file = pic_file, width=300, height=500)
plot(sort(sums, decreasing=T), pch=19, col="green", main="OTU abundance distribution", xlab="Index", ylab="read count", type="p")
points(1:length(sums), rep(10000, as.integer(length(sums))), col=c("black"), type="l", lwd=2)
dev.off()

# Normalize by sums
normalized_clean_otus <- matrix(data=NA, nrow=NROW(clean_otus), ncol=NCOL(clean_otus))
for (i in 1:NCOL(clean_otus)) {
  normalized_clean_otus[,i] <- clean_otus[,i] / sums[i]
}

# Put them back into the original dataframe (with col/row headers)  
normalized_otu_table <- as.matrix(otu_table_removed)
normalized_otu_table[1:NROW(normalized_otu_table),2:(NCOL(normalized_otu_table)-1)] <- as.numeric(normalized_clean_otus)

# Sanity check to make sure they still sum up to 1
sum(as.numeric(normalized_otu_table[1:NROW(normalized_otu_table),2]))

# Remove any OTUs with less than 10000 reads
normalized_otu_table_10000 <- normalized_otu_table[,-(which(sums<10000))]

# Remove samples that have less than 2 replicates and take the average of the samples that have 2 or more
sample_names <- names(table(colnames(normalized_otu_table_10000))[table(colnames(normalized_otu_table_10000)) > 1])
OTU_names <- as.character(normalized_otu_table_10000[,1])
num_row <- NROW(normalized_otu_table_10000)
OTU_Averages <- matrix(data=NA, ncol=length(sample_names), nrow=num_row, dimnames=list(OTU_names, sample_names))
for (i in 1:length(sample_names)) { 
  pos <- which(colnames(normalized_otu_table_10000) %in% sample_names[i])
  accession_matrix <- matrix(data=as.numeric(normalized_otu_table_10000[,pos]), nrow=num_row, ncol=length(pos))
  averages <- apply(accession_matrix, 1, sum) / length(pos)
  OTU_Averages[,i] <- averages
}
heatmap(OTU_Averages)

# Split the taxonomy
tax_split <- strsplit(as.character(tax1$taxonomy), "; ")
parts_matrix <- matrix(unlist(tax_split), ncol=length(tax_split[[1]]), byrow = TRUE)
split_tax <- cbind(parts_matrix, as.character(tax1$id))
split_tax <- as.data.frame(split_tax, stringssAsFactors=F)
colnames(split_tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")

# Sum read counts by group per taxon level
Averages <- c()
Averages[[1]] <- aggregate(OTU_Averages, by=list(split_tax[,1]), FUN=sum)
Averages[[2]] <- aggregate(OTU_Averages, by=list(split_tax[,2]), FUN=sum)
Averages[[3]] <- aggregate(OTU_Averages, by=list(split_tax[,3]), FUN=sum)
Averages[[4]] <- aggregate(OTU_Averages, by=list(split_tax[,4]), FUN=sum)
Averages[[5]] <- aggregate(OTU_Averages, by=list(split_tax[,5]), FUN=sum)
Averages[[6]] <- aggregate(OTU_Averages, by=list(split_tax[,6]), FUN=sum)
Averages[[7]] <- aggregate(OTU_Averages, by=list(split_tax[,7]), FUN=sum)

# Add the row names as values
Averages[[8]] <- cbind(rownames(OTU_Averages), OTU_Averages)

# Restructure the genotype file
marker_table <- read.table(file.path(inp_folder, "markers_new_map.csv"), sep = ",", comment.char="", header=F, check.names=F, stringsAsFactors=F)
marker_table[1,1] <- "id"
marker_table[2,1] <- ""
marker_table[3,1] <- ""
marker_table <- marker_table[-4,]
marker_table[marker_table == "BB"] <- "B"
marker_table[marker_table == "AA"] <- "A"

# Performs the analysis for each taxonomy level
for (taxon_level in 6:6) {
  ave_tax <- Averages[[taxon_level]]
  for (i in 2:NCOL(ave_tax)) {
    ave_tax[,i] <- as.numeric(ave_tax[,i])
  }
  
  # Add the column names as values
  A_updated <- rbind(colnames(ave_tax), ave_tax)
  
  # Transpose the table
  t_Averages <- t(A_updated)
  t_Averages[1,1] <- "id"
  rownames(t_Averages)[[1]] <- "id"
  
  # Write relative read counts of the current taxonomic level to file
  count_file <- file.path(inp_folder, sprintf("%s_relative_readcounts_taxlevel%s.txt", prefix, taxon_level))
  write.table(t_Averages, file=count_file, sep=",", quote=F, col.names=F, row.names=F)
  
  # Remove samples that are missing in either the genotype or the phenotype
  phenotype <- t_Averages[which(rownames(t_Averages) %in% marker_table[,1]),]
  genotype <- marker_table[which(marker_table[,1] %in% c("", t_Averages[,1])),]
  
  # Keep the top 500 taxons with the largest SD (if there are more than 500 taxons)
  OTU_SD <- apply(phenotype[2:NROW(phenotype),2:NCOL(phenotype)], 2, sd)
  if (length(OTU_SD) > 500) {
    top500_SD_cut_off <- as.numeric(sort(OTU_SD, decreasing=TRUE)[500])
    top500_SD_positions <- which(as.numeric(OTU_SD) >= top500_SD_cut_off)
    top500_SD_positions <- c(1, top500_SD_positions+1)
    phenotype <- phenotype[,top500_SD_positions]
    # Keep the taxonomy of the top 500
    split_tax <- split_tax[which(split_tax[[colnames(split_tax)[[taxon_level]]]] %in% colnames(top500_SD_OTU_Averages)),]
  }
  
  # Write the tables to the output files
  tax_file <- file.path(inp_folder, sprintf("%s_taxonomy_taxlevel%s.txt", prefix, taxon_level))
  write.table(tax1, file=tax_file, sep=",", quote=F, col.names=F, row.names=F)
  otu_file <- file.path(inp_folder, sprintf("%s_phenotype_taxlevel%s.txt", prefix, taxon_level))
  write.table(phenotype, file=otu_file, sep=",", quote=F, col.names=F, row.names=F)
  marker_file <- file.path(inp_folder, sprintf("%s_genotype_taxlevel%s.txt", prefix, taxon_level))
  write.table(genotype, file=marker_file, sep=",", quote=F, col.names=F, row.names=F)
  
  if (NCOL(phenotype) < 3) {
    next
  }
  
  # Read files and transform into a readable format for R/qtl2
  genfile <- file.path(inp_folder, sprintf("%s_genotype_taxlevel%s.txt", prefix, taxon_level))
  phefile <- file.path(inp_folder, sprintf("%s_phenotype_taxlevel%s.txt", prefix, taxon_level))
  in_cross <- qtl::read.cross("csvs", genfile = genfile, phefile = phefile, crosstype = "riself", sep=",")
  mycross <- (qtl2::convert2cross2(in_cross))
  # Creates markers at the given interval
  #map <- insert_pseudomarkers(mycross$gmap, step=1)
  # Calculates the probability of a genotype at each marker
  pr <- calc_genoprob(mycross, mycross$gmap, error_prob=0.0001)
  # Permutation of genome scan using Haley-Knott regression
  out_perm <- scan1perm(pr, mycross$pheno, n_perm=10000)

  # Contains highest LOD score for each phenotype across the entire genome
  perm_thres <- summary(out_perm)
  if ("id" %in% colnames(perm_thres)) {
    perm_thres <- perm_thres[,-1]
  }
  perm_thres[perm_thres == -Inf] <- 0
  # out_perm contains the maximum LOD score for each OTU for each permutation.
  # Genome scan using Haley-Knott regression
  out <- scan1(pr, mycross$pheno)
  out[is.nan(out)] <- 0
  
  if ("id" %in% colnames(out)) {
    out <- out[,-1]
  }
  
  # Restructure the map matrix
  map_updated <- c("Chromosome", "cM")
  for (i in 1:NROW(mycross$gmap)) {
    for (j in 1:NROW(mycross$gmap[[i]])) {
      map_updated <- rbind(map_updated, cbind(i, mycross$gmap[[i]][[j]]))
    }
  }
  
  # Add row and column names and mapping locations to LOD matrix
  out_update <- as.data.frame(out)
  out_updated <- cbind(rownames(out_update), as.data.frame(out_update, stringsAsFactors=F), stringsAsFactors=F)
  out_updatede <- rbind(c("SNP.name", colnames(out_update), stringsAsFactors=F), out_updated, stringsAsFactors=F)
  out_updateded <- cbind(map_updated, out_updatede, stringsAsFactors=F)
  
  # Write results to file
  out_file <- file.path(inp_folder, sprintf("%s_lods_taxlevel%s.txt", prefix, taxon_level))
  write.table(out_updateded, file=out_file, sep="\t", quote=F, col.names=F, row.names=F)
  
  # If there is only 1 column with LOD scores, there can't be peaks
  if (NCOL(out) < 2) {
    next
  }
  
  # Find significant LOD peaks
  threshold <- quantile(perm_thres, c(.95))
  print(threshold)
  thres_file <- file.path(inp_folder, sprintf("%s_threshold_taxlevel%s.txt", prefix, taxon_level))
  write.table(threshold, file=thres_file, sep="\t", quote=F, col.names=F, row.names=F)
  peaks <- as.data.frame(find_peaks(out, mycross$gmap, threshold=threshold, peakdrop=0), stringsAsFactors=F)
  print(peaks)
  print(taxon_level)
  # Restructure the peaks matrix
  peaks_updated <- cbind(peaks, peaks[,2], stringsAsFactors=F)
  colnames(peaks_updated)[[2]] <- "id"
  peaks_updatede <- as.data.frame(merge(peaks_updated, split_tax, by.x="id", by.y=colnames(split_tax)[[taxon_level]], sort=F), stringsAsFactors=F)
  index_num <- taxon_level + 5
  peaks_updateded <- peaks_updatede[,c(3:index_num)]
  colnames(peaks_updateded)[[4]] <- colnames(split_tax)[[taxon_level]]
  taxon_level_colnames <- colnames(split_tax)[1:taxon_level]
  peaks_updateded <- peaks_updateded[c(colnames(peaks_updateded)[1:3], taxon_level_colnames)]
  peaks_updatedede <- rbind(as.character(colnames(peaks_updateded)), as.data.frame(peaks_updateded, stringsAsFactors=F), stringsAsFactors=F)
  for (i in 1:NCOL(peaks_updatedede)) {
    peaks_updatedede[,i] <- as.character(peaks_updatedede[,i])
  }
  peaks_updatedede[1,] <- c("Chromosome", "cM", "LOD_score", taxon_level_colnames)
  peaks_updatedede <- unique(peaks_updatedede)
  # Write results to files
  sig_file <- file.path(inp_folder, sprintf("%s_peaks_taxlevel%s.txt", prefix, taxon_level))
  write.table(peaks_updatedede, file=sig_file, sep="\t", quote=F, col.names=F, row.names=F)
  
  # If there are no peaks, there are no proteins to be found
  if (length(peaks[[2]]) < 1) {
    next
  }
  
  # File paths
  lod_file <- file.path(inp_folder, sprintf("%s_lods_taxlevel%s.txt", prefix, taxon_level))
  peak_file <- file.path(inp_folder, sprintf("%s_peaks_taxlevel%s.txt", prefix, taxon_level))
  
  # Read files
  lods <- read.csv(file = lod_file, sep = '\t', header = TRUE)
  peaks <- read.csv(file = peak_file, sep = '\t', header = TRUE)
  
  ymax <- max(lods[,5:NCOL(lods)])
  
  custom_palette <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
  custom_palette[custom_palette == "#999999"] <- "black"
  custom_palette[custom_palette == "#B3B3B3"] <- "cyan"
  palette(custom_palette)
  colours <- c()
  
  for (num in 1:12) {
    sub_peaks <- subset(peaks, Chromosome == num)
    # if (length(sub_peaks$cM) == 0) {
    #   next
    # }
    pic_file <- file.path(inp_folder, sprintf("%s_lods_taxlevel%s_chromosome_%s.jpeg", prefix, taxon_level, num))
    jpeg(file = pic_file, width=1200, height=600)
    par(mar=c(4, 4, 4, 8), xpd=TRUE)
    sub_lods <- subset(lods, Chromosome == num)
    plot(sub_lods$cM, sub_lods[[4]], pch=19, col="gray", main=sprintf("LOD curves chromosome %s", num), ylim=c(0, ymax), xlab="cM", ylab="LOD score", type="l", cex=5)
    for (i in 5:NCOL(sub_lods)) {
      points(sub_lods$cM, sub_lods[[i]], pch=19, col="gray", type="l")
    }
    names.use <- names(sub_lods)[(names(sub_lods) %in% sub_peaks[[colnames(split_tax)[[taxon_level]]]])]
    fctr <- factor(names.use)
    cc <- c(fctr)
    cc[cc == "8"] <- "orange"
    if (length(names.use) > 0) {
      names.use <- c(c("SNP.name", "Chromosome", "cM"), names.use)
      otu_set <- sub_lods[, names.use]
      for (i in 4:NCOL(otu_set)) {
        col_num <- i-3
        points(otu_set[[3]], otu_set[[i]], pch=19, col=cc[[col_num]], type="l")
      }
    }
    fctr2 <- factor(sub_peaks[[colnames(split_tax)[[taxon_level]]]])
    cc2 <- c(fctr2)
    cc2[cc2 == "8"] <- "orange"
    points(sub_peaks$cM, sub_peaks$LOD_score, pch=19, col=cc2)
    points(sub_lods$cM, rep(as.numeric(threshold), as.integer(length(sub_lods$cM))), col=c("black"), type="l", lwd=2)
    colours <- c(colours, cc)
    if (length(fctr) > 0) {
      legend("topright", inset=c(0,0), legend=unique(fctr), pch=19, col=unique(cc), title="Name")
    } else {
      legend("topright", inset=c(0,0), legend=c("None"), pch=19, col=unique(cc), title="Name")
    }
    dev.off()
  }
  
  otu_names <- c()
  for (i in 1:NROW(peaks[[2]])) {
    otu_names <- c(otu_names, paste0(paste(replicate(90, " "), collapse = ""), peaks[[2]][[i]]))
  }
  n_and_t <- c(as.character(peaks$taxonomy))
  coloration <- c(colours, c(rep(0, length(n_and_t))))
  n_and_t <- c(as.character(otu_names), n_and_t)
  pic_file <- file.path(inp_folder, sprintf("%s_legend_taxlevel%s.jpeg", prefix, taxon_level))
  jpeg(file = pic_file, width=1200, height=600)
  par(mar=c(1, 1, 1, 1), xpd=TRUE)
  plot(0,type='n',axes=FALSE,ann=FALSE)
  legend("topright", col=coloration, ncol=2, legend=n_and_t, pch=19, x.intersp=-22, title="              OTU name           Taxonomy", title.adj=0)
  dev.off()
  
  # Read protein file and discard proteins that aren't assigned to a chromosome
  protein_file <- file.path(inp_folder, "ProteinTable7_211970.txt")
  proteins <- read.csv(file = protein_file, sep = '\t', header = TRUE)
  proteins <- subset(proteins, X.Replicon.Name != "Un")
  
  # Get basepair positions of markers
  bp <- c()
  for (i in unique(lods$Chromosome)) {
    sub_lods <- subset(lods, Chromosome == i)
    for (num in 1:NROW(sub_lods)) {
      # First marker of a chromosome is always on 0 cM -> 0 bp is used
      if (num == 1) {
        bp <- c(bp, 0)
        next
      }
      #
      lod <- sub_lods[num,]
      # Real markers have their basepair position in their name
      # Those positions are used to infer the position of the other markers
      if (grepl("-", lod$SNP.name)) {
        name <- strsplit(as.character(lod$SNP.name), "-")
        first_part <- as.character(name[[1]][[1]])
        # Except for 3 markers with 'seq' in their name
        # Their positions are inferred in the same manner as the pseudo-markers
        if (first_part != "seq") {
          bp <- c(bp, as.integer(first_part))
          next
        }
      }
      # Find closest real marker on the left (lower cM/bp)
      low_sub <- subset(sub_lods, cM < as.numeric(lod$cM))
      semi_sub <- subset(low_sub, grepl("-", low_sub$SNP.name))
      real_sub <- subset(semi_sub, !grepl("seq", semi_sub$SNP.name))
      low_idx <- which.max(real_sub$cM)
      lower <- real_sub[low_idx,]
      # Find closest real marker on the right (higher cM/bp)
      high_sub <- subset(sub_lods, cM > as.numeric(lod$cM))
      semi_sub <- subset(high_sub, grepl("-", high_sub$SNP.name))
      real_sub <- subset(semi_sub, !grepl("seq", semi_sub$SNP.name))
      high_idx <- which.min(real_sub$cM)
      higher <- real_sub[high_idx,]
      # Determine position of marker in basepairs
      low_name <- strsplit(as.character(lower$SNP.name), "-")
      low_bp <- as.integer(low_name[[1]][[1]])
      high_name <- strsplit(as.character(higher$SNP.name), "-")
      high_bp <- as.integer(high_name[[1]][[1]])
      distance_bp <- high_bp - low_bp
      distance_cm <- as.numeric(higher$cM) - as.numeric(lower$cM)
      multiplication_factor <- distance_bp / distance_cm
      diff_cm <- lod$cM - lower$cM
      lod_bp <- (multiplication_factor * diff_cm) + low_bp
      # Add position to list
      bp <- c(bp, lod_bp)
    }
  }
  # Add basepair position list to the rest of the data
  lods <- cbind(lods, bp)
  
  # Get search space for each significant marker (bound by neighbouring markers)
  boundaries <- c()
  for (num in unique(lods$Chromosome)) {
    sub_lods <- subset(lods, Chromosome == num)
    for (i in 1:NROW(sub_lods)) {
      lod <- sub_lods[i,]
      for (n in 1:NROW(peaks)) {
        peak <- peaks[n,]
        if (lod$Chromosome == peak[[1]] & lod$cM == peak[[2]]) {
          for (j in i:1) {
            if (sub_lods[j,][[as.character(peak[[colnames(split_tax)[[taxon_level]]]])]] < threshold) {
              LL <- as.integer(sub_lods[j,]$bp)
              break
            }
          }
          for (j in i:NROW(sub_lods)) {
            if (sub_lods[j,][[as.character(peak[[colnames(split_tax)[[taxon_level]]]])]] < threshold) {
              UL <- as.integer(sub_lods[j,]$bp)
              break
            }
          }
          #LL <- as.integer(lods[i-1,]$bp)
          #UL <- as.integer(lods[i+1,]$bp)
          boundaries <- rbind(boundaries, c(LL, UL))
        }
      }
    }
  }
  
  # Attach boundaries to peaks (a.k.a. significant markers)
  peaks <- cbind(peaks, boundaries)
  
  # Get proteins within search space
  possible_prots <- c()
  for (i in 1:12) {
    sub_prots <- subset(proteins, X.Replicon.Name == i)
    sub_peaks <- subset(peaks, Chromosome == i)
    if (NROW(sub_peaks) == 0) {
      next
    }
    for (n in 1:NROW(sub_peaks)) {
      peak <- sub_peaks[n,]
      # Protein start position lies within boundaries
      start_set <- subset(sub_prots, Start > peak$"1" & Start < peak$"2")
      # Protein stop position lies within boundaries
      stop_set <- subset(sub_prots, Stop > peak$"1" & Stop < peak$"2")
      # Boundaries lie within protein start and stop
      within_set <- subset(sub_prots, Start > peak$"1" & Stop < peak$"2")
      # Discard doubles
      unified <- merge(start_set, stop_set, by="Protein.product")
      if (NROW(within_set) != 0) {
        unified <- merge(unified, within_set, by="Protein.product")
      }
      # Bind proteins to peaks
      if (NROW(unified) == 0) {
        next
      } else if (NROW(unified) == 1) {
        combo <- c(peak, unified[1,])
        possible_prots <- rbind(possible_prots, combo)
      } else {
        for (i in 1:NROW(unified)) {
          combo <- c(peak, unified[i,])
          possible_prots <- rbind(possible_prots, combo)
        }
      }
    }
  }
  
  # Write peak-protein combinations to a file
  peak_prots_file <- file.path(inp_folder, sprintf("%s_peak_proteins_taxlevel%s.txt", prefix, taxon_level))
  if (file.exists(peak_prots_file)) file.remove(peak_prots_file)
  write.table(possible_prots[1,], file=peak_prots_file, sep="\t", append=TRUE, quote=FALSE)
  for (i in 2:NROW(possible_prots)){
    write.table(possible_prots[i,], file=peak_prots_file, sep="\t", append=TRUE, quote=FALSE, col.names=FALSE)
  }
  
  # The part below is for exporting the marker-OTU combinations to Cytoscape
  
  library(devtools)
  #install_github('cytoscape/RCy3')
  library(RCy3)
  
  
  nodes <- matrix(nrow=0, ncol=4)
  colnames(nodes) <- c("id", "colour_idx", "chromosome", "cM")
  edges <- matrix(nrow=0, ncol=2)
  colnames(edges) <- c("source","target")
  
  for (i in 1:NROW(possible_prots)) {
    poss_prot <- possible_prots[i,]
    nodes <- rbind(nodes, c(as.character(poss_prot[[colnames(split_tax)[[taxon_level]]]]), "1", "", ""))
    nodes <- rbind(nodes, c(as.character(poss_prot$Protein.name.x), "2", as.character(poss_prot$Chromosome), as.character(poss_prot$cM)))
    edges <- rbind(edges, c(as.character(poss_prot[[colnames(split_tax)[[taxon_level]]]]), as.character(poss_prot$Protein.name.x)))
  }
  nodes <- unique(nodes)
  edges <- unique(edges)
  
  createNetworkFromDataFrames(data.frame(nodes, stringsAsFactors = F),
                              data.frame(edges, stringsAsFactors = F),
                              title=sprintf("%s_taxlevel%s.txt", prefix, taxon_level),
                              collection="Network Traits and Genomes")
}


time.taken <- Sys.time() - start.time
print("Result files written and script completed.")
time.taken
