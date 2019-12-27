#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 4/4/17
#
# eQTL_network.R
#
# Goal is to build a regulatory network based on eQTL
#
# Usage:
# eQTL_summary.R EQTL=fdr.trans.eqtls.txt [Options]
#  EQTL=    Specifies the main input file of EQTLs from eQTL_perm_fdr.R
#  OUTPUT=  Path and file stub for all output files (defaults to EQTL minus the .trans.eqtls.txt suffix)
#  CORES=   Number of threads to use for parallelization
#
# The output format is in the format:
# GENEX GENEY [Details]
# Meaning that GeneX is inferred to regulate GeneY, with details showing the SNPs involved, etc.
#
# NOTE: The output format MAY CHANGE in the future, depending on what tools are used downstream of this step!
#
# Currently this only uses GeneX->GeneY connections where the same SNP is a cis eQTL for GeneX and a trans eQTL for GeneY
# TO DO: Add an option to also include cases where a trans eQTL for GeneY is merely proximal to GeneX?
#

# TO DO: A major issue is cycles within the regulatory network, which complicates the process of determining number of ancestors/dependents in the graph
# I have 2 ideas for solutions:
#
# Solution 1: [NOT USED]
# Compute total upstream/downstream in parallel, as follows:
# For each gene X:
#   Compute the set of immediate upstream genes, all u in U s.t. u->X
#   Record |U| as DIRECT.UPSTREAM
#   Compute the set of immediate downstream genes, all d in D s.t. X->d AND d not in U
#   Record |D| as DIRECT.DOWNSTREAM
#   Repeat until |U| and |D| both stop changing:
#     Add to U all genes v s.t. v->u for some u in U AND v not in D
#     Add to D all genes c s.t. d->c for some d in D AND c not in U
#   Record |U| as TOTAL.UPSTREAM
#   Record |D| as TOTAL.DOWNSTREAM
#
# Intuitively, this should enforce that each other gene Y!=X in the graph is EITHER upstream or downstream of X, but NEVER both
# And in the presence of loops, each gene is put in the category for which it has the shortest directed path to X
# In the presence of a tie, that is X->...->Y and Y->...->X are the same length, then Y is considered upstream
#
# The advantage of this approach is we don't have to drop any edges, we can preserve the original network with cycles
# The disadvantage is that we are still treating all edges as equally good, which may not be the case
#
#
# Solution 2: [IMPLEMENTED BELOW]
# First enforce that the eQTL network be acyclic before proceeding to the computation of total upstream/downstream nodes for each gene
# An efficient way to do this for each disjoint subnetwork, is the following:
#   Rank all the edges in the subnetwork from strongest to weakest (easiest way is to take the min p-value for any SNP, and rank edges from lowest to highest of these p-values)
#     In any case, this gives a ranking of all edges, E_i for i in 1:N
#   Start a new DAG subnetwork with E_1 only (that is, the most significant edge), connecting X_1->Y_1, denote this set of edges E'
#   Set depth X_1=0, depth Y_1=1, depth of all other genes = NA (not yet in E')
#   For each i in 2:N, consider E_i connecting new X_i->Y_i:
#     If depth of X_i and Y_i are both defined, and depth Y_i < depth X_i:
#       Remove E_i from the network (or put in set E_0)
#     Else add E_i to E', then recompute the depth of all genes in E' as follows:
#       Identify all roots R, consisting of each node r for which there is no edge directed into r
#       Set depth of all roots R to 0, set index d = 0
#       Put all other nodes in E' in the undetermined set U
#       Repeat until U is empty:
#         Identify all direct children C of nodes with depth d
#         For all nodes C, set depth = d+1
#         Remove all nodes C from U
#         Increment d
#       (After this loop, all nodes currently in E' are given their minimal depth from any root node)
#  (After this loop, all edges in E have either been added to E' or E_0. E' is a true DAG containing the highest ranking acyclic edges of E, |E_0| = number of low ranking edges that needed to be dropped)
#  
   


# --- INITIALIZATION --- #

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")

# For parallelization
library(foreach)
library(doMC)

options(stringsAsFactors=F)

usageStr="USAGE:\neQTL_summary.R EQTL=fdr.trans.results.txt MAP=snp.gene.map"
my.args <- commandArgs(trailingOnly=T)
# TEMP TESTING:
# my.args <- c("EQTL=known_all_novel_genes/plink/ExprF.0.05.fdr.trans.eqtls.txt")

if(length(my.args) == 0) {
  stop("Requires at least one argument. ", usageStr)
} else {
  cat("Running with params:",my.args,"\n")
}
# Parse param names and values
argSplit=strsplit(my.args, split="=", fixed=T)
my.args <- lapply(argSplit, "[", 2)
names(my.args) <- unlist(lapply(argSplit, "[", 1))

# Fill in default value for missing params
setDefault <- function(param, default) {
  if(!(param %in% names(my.args))) {
    cat("Setting ",param,"=",default," by default\n", sep="")
    my.args[[param]] <<- default
  }
}

# EQTL=  The file containing line mean expression values to filter and prep for downstream analysis
# Make sure EQTL param exists
if(!("EQTL" %in% names(my.args))) {
  stop("Missing EQTL parameter. ", usageStr)
}

default.outstub <- sub("[.]txt$", "", my.args$EQTL)
default.outstub <- sub("[.]eqtls$", "", default.outstub)
default.outstub <- sub("[.]trans$", "", default.outstub)
setDefault("OUTPUT",default.outstub)

outDir <- dirname(my.args$OUTPUT)
if(!file.exists(outDir)) {
  dir.create(path=outDir, recursive=T, showWarnings = F)
}

# CORES= Number of CPUs to use for multi-threaded tasks
setDefault("CORES",as.integer(NA))
if(is.na(my.args$CORES)) {
  # If NA, see if environment variable SLURM_JOB_CPUS_PER_NODE exists
  slurm.cpus <- Sys.getenv("SLURM_JOB_CPUS_PER_NODE")
  if(slurm.cpus != "") {
    slurm.cpus <- as.integer(slurm.cpus)
    if(slurm.cpus > 0) {
      my.args$CORES <- slurm.cpus
      cat("Inferring CORES =", slurm.cpus, "based on SLURM_JOB_CPUS_PER_NODE.\n")
    }
  }
}
# If STILL NA, use 1 by default
if(is.na(my.args$CORES)) {
  my.args$CORES <- 1
  cat("Setting CORES = 1 by default.\n")
}
# If CORES = ALL, use detectCores()-1
if(is.character(my.args$CORE) & (my.args$CORE == "ALL")) {
  my.args$CORES <- as.integer(detectCores()-1)
  cat("Setting CORES =", my.args$CORES, "\n")
}
# Otherwise, attempt to convert to integer
if(!is.integer(my.args$CORE)) {
  my.args$CORE <- as.integer(my.args$CORE)
  cat("Setting CORES =", my.args$CORE, "\n")
}

# Set up parallel backend with desired number of CPUs
cat("Using", my.args$CORES, "CPUs for parallel tasks.\n")
registerDoMC(my.args$CORES)
cat("\n")


# --- Load Input Files --- #

# Load the EQTL File
eqtl.table <- read.table(my.args$EQTL, header=T, sep="\t")
cat("Loaded trans eQTL results for", nrow(eqtl.table), "trans eQTLs from", my.args$EQTL, "\n")

# Drop all features with no CIS flags
drop.features <- !grepl("CIS", eqtl.table$SNP.GENE.REASON)
if(any(drop.features)) {
  cat("Dropping", sum(drop.features), "trans eQTLs with no cis gene targets.\n")
  eqtl.table <- eqtl.table[!drop.features,]
  cat(nrow(eqtl.table), "trans eQTLs remain for network analysis.\n")
}
cat("\n")

# Testing: reduce to just the first 100 trans eQTLs
# eqtl.table <- eqtl.table[1:1000,]

# Expand the table so that there is a separate row for every Gene,SNP,Gene triplet (split out the cases with multiple potential upstream genes)
# TO DO: Could easily parallelize this
exp.eqtl.table <- foreach(i=1:nrow(eqtl.table), .combine=rbind) %dopar% {
  exp.n <- length(unlist(strsplit(eqtl.table[i,"SNP.GENE"], split=",", fixed=T)))
  if(exp.n==1) {
    return(eqtl.table[i,])
  } else {
    exp.rows <- data.frame(row.names=1:exp.n)
    for(j in colnames(eqtl.table)) {
      if(grepl("[,]",eqtl.table[i,j])) {
        exp.rows[,j] <- unlist(strsplit(eqtl.table[i,j], split=",", fixed=T))
      } else {
        exp.rows[,j] <- rep(eqtl.table[i,j], times=exp.n)
      }
    }
    row.names(exp.rows) <- NULL
    # Make sure we're only keeping cis eQTL-supported connections
    exp.rows <- exp.rows[exp.rows$SNP.GENE.REASON=="CIS",]
    # TO DO: Could try to choose best case here, e.g. prioritizing 5PRIME > EXON > INTRON > 3PRIME
    return(exp.rows)
  }
}

cat("Expanded", nrow(eqtl.table),"trans eQTL SNPs to", nrow(exp.eqtl.table), "Gene,SNP,Gene trios.\n\n")
eqtl.table <- exp.eqtl.table
rm(exp.eqtl.table)

# Now get list of unique SNP.GENE,GENE pairs
gene.pairs <- unique(paste(eqtl.table$SNP.GENE,eqtl.table$GENE,sep=","))
cat("Summarizing all SNPs for",length(gene.pairs),"gene->gene connections.\n")
gene.pairs <- strsplit(gene.pairs, split=",", fixed=T)
names(gene.pairs) <- unlist(lapply(gene.pairs, paste, collapse=","))

# Now for each pair, build a table of all SNPs
cat("Extracting SNPs for each gene->gene connections.\n")
gene.pair.snps <- foreach(gp=gene.pairs) %dopar% {
  eqtl.table[(eqtl.table$SNP.GENE==gp[1])&(eqtl.table$GENE==gp[2]),]
}
names(gene.pair.snps) <- names(gene.pairs)

# Build the gene.pair table with some summary columns of the SNPs linking them together
cat("Building table for all gene->gene connections.\n")
gene.pair.table <- foreach(gp=names(gene.pairs), .combine=rbind) %dopar% {
  gp.snps <- gene.pair.snps[[gp]]
  # Make sure SNPs sorted by significance
  gp.snps <- gp.snps[order(gp.snps$PVAL, decreasing=F),]
  data.frame(REG.GENE=gene.pairs[[gp]][1], TARGET.GENE=gene.pairs[[gp]][2], SNP=paste(gp.snps$SNP, collapse=","), SNP.PVAL=paste(gp.snps$PVAL, collapse=","), SNP.REG.DIST=paste(gp.snps$SNP.GENE.DIST, collapse=","), SNP.REG.PART=paste(gp.snps$SNP.GENE.PART, collapse=","))
}

# TO DO: Could load coordinates of both genes from gene info table
# TO DO: Could load symbols of genes here

# Identify disparate (unconnected) sub-networks
# e.g. take the first gene and add all genes connected in either direction,
# repeat until there are no new edges to add,
# then repeat the whole process on remaining genes
gene.pair.table[,"NETWORK"] <- as.integer(NA)
unlabeled.genes <- unique(c(gene.pair.table$REG.GENE, gene.pair.table$TARGET.GENE))
cat("Splitting", length(unlabeled.genes), "genes into distinct subnetworks.\n")
next.net <- 1
while(length(unlabeled.genes) > 0) {
  # Seed with next gene not assigned to a network
  network.genes <- unlabeled.genes[1]
  network.edges <- rep(F, times=nrow(gene.pair.table))
  # Loop until network.edges is stable
  repeat {
    # Identify all edges encompassing current gene list
    exp.network.edges <- (gene.pair.table$REG.GENE %in% network.genes) | (gene.pair.table$TARGET.GENE %in% network.genes)
    # These should all be unassigned
    stopifnot(all(is.na(gene.pair.table$NETWORK[exp.network.edges])))
    if(all(network.edges == exp.network.edges)) {
      break
    } else {
      network.edges <- exp.network.edges
      network.genes <- unique(c(gene.pair.table[network.edges,c("REG.GENE","TARGET.GENE")],recursive=T))
    }
  }
  gene.pair.table[network.edges,"NETWORK"] <- next.net
  next.net <- next.net+1
  # Remove all network genes from unlabeled.genes
  unlabeled.genes <- setdiff(unlabeled.genes, network.genes)
}
# How many networks?
cat("Identified",length(unique(gene.pair.table$NETWORK)),"distinct sub-networks.\n\n")
# There should be no unassigned edges
stopifnot(sum(is.na(gene.pair.table$NETWORK))==0)

# To Pick up here:
# gene.pair.table <- read.table(paste0(my.args$OUTPUT, ".cis.trans.network.txt"), header=T, sep="\t")

# Technically, this loop could be parallelized, but there's no real need
# The main network will take the longest, and everything else should be instanteous
# Unfortunately, I don't see a way to parallelize the inner loops
gene.pair.filtered <- foreach(net=1:max(gene.pair.table$NETWORK), .combine='rbind') %do% {
  gene.pair.subnet <- gene.pair.table[gene.pair.table$NETWORK==net,]
  
  # First enforce that the eQTL network be acyclic before proceeding to the computation of total upstream/downstream nodes for each gene
  # An efficient way to do this for each disjoint subnetwork, is the following:
  #   Rank all the edges in the subnetwork from strongest to weakest (easiest way is to take the min p-value for any SNP, and rank edges from lowest to highest of these p-values)
  #     In any case, this gives a ranking of all edges, E_i for i in 1:N
  
  # Unpack the SNP P-vals, keep the best (min)
  gene.pair.subnet[,"SNP.BEST.PVAL"] <- unlist(lapply(lapply(strsplit(gene.pair.subnet$SNP.PVAL, split=","), as.numeric), min))
  # Rank by min P-val
  gene.pair.subnet <- gene.pair.subnet[order(gene.pair.subnet$SNP.BEST.PVAL, decreasing=F),]
  
  # Start a new DAG subnetwork E', initialize with no edges
  # Storing this as just a boolean vector indicating if each edge from E belongs in E'
  # All initialized as F, to make it easy to get the current E' as we're building
  # e.g. E' = E[subnet.dag,]
  subnet.dag <- rep(F, times=nrow(gene.pair.subnet))
  # Set depth of all edges in E as -1 (not yet in E')
  subnet.genes <- unique(c(gene.pair.subnet$REG.GENE, gene.pair.subnet$TARGET.GENE))
  subnet.gene.depth <- rep(as.integer(-1), times=length(subnet.genes))
  names(subnet.gene.depth) <- subnet.genes
  # For each i in 1:N, consider E_i connecting new X_i->Y_i:
  # This loop CANNOT be parallelized!
  for(i in 1:nrow(gene.pair.subnet)) {
    xi <- gene.pair.subnet[i,"REG.GENE"]
    yi <- gene.pair.subnet[i,"TARGET.GENE"]
    if(all(subnet.gene.depth[c(xi,yi)] != -1) & (subnet.gene.depth[yi] < subnet.gene.depth[xi])) {
      # If depth of X_i and Y_i are both defined, and depth Y_i < depth X_i:
      # Remove E_i from the network
      subnet.dag[i] <- F
    } else {
      # Else add E_i to E', then recompute the depth of all genes in E'
      subnet.dag[i] <- T
      cur.subnet <- gene.pair.subnet[subnet.dag,]
      # Recompute depth of all nodes in E':
      # Identify all roots R, consisting of each node r for which there is no edge directed into r
      cur.roots <- setdiff(cur.subnet$REG.GENE, cur.subnet$TARGET.GENE)
      # Set depth of all roots R to 0
      subnet.gene.depth[cur.roots] <- 0
      # set current depth = 0
      cur.depth <- 0
      # Put all other nodes in E' in the undetermined set U
      undeterm.genes <- setdiff(unique(c(cur.subnet$REG.GENE, cur.subnet$TARGET.GENE)), cur.roots)
      # Repeat until U is empty:
      while(length(undeterm.genes) > 0) {
        # Identify all direct children C of nodes with depth d
        cur.parents <- names(subnet.gene.depth)[subnet.gene.depth==cur.depth]
        cur.children <- unique(cur.subnet[cur.subnet$REG.GENE %in% cur.parents,"TARGET.GENE"])
        # For all nodes C, set depth = d+1
        subnet.gene.depth[cur.children] <- cur.depth + 1
        # Remove all nodes C from U
        undeterm.genes <- setdiff(undeterm.genes, cur.children)
        # Increment d
        cur.depth <- cur.depth + 1
      }
      # (After this loop, all nodes currently in E' are given their minimal depth from any root node)
    }
  }
  
  # All edges in E have either been added to E' (T) or E_0 (F)
  # E' is a true DAG containing the highest ranking acyclic edges of E
  # |E_0| = number of low ranking edges that needed to be dropped)
  cat("Dropping", sum(!subnet.dag), "edges from subnetwork", net, "to break loops.\n")
  cat("DAG subnetwork", net, "has", sum(subnet.dag), "edges remaining.\n\n")
  return(gene.pair.subnet[subnet.dag,])
}

# Replace full gene pair table with the filtered one:
gene.pair.table <- gene.pair.filtered[,colnames(gene.pair.table)]
rm(gene.pair.filtered)

# Re-number the networks by size (# of edges)
network.edge.sz <- table(gene.pair.table$NETWORK)
network.order <- 1:max(gene.pair.table$NETWORK)
names(network.order) <- order(network.edge.sz, decreasing=T)
gene.pair.table$NETWORK <- network.order[as.character(gene.pair.table$NETWORK)]
gene.pair.table <- gene.pair.table[order(gene.pair.table$NETWORK, decreasing=F),]
network.edge.sz <- table(gene.pair.table$NETWORK)

# Report size of each network (edges and genes)
cat("Range of network sizes (by edge #):", range(network.edge.sz), "\n")
cat("Mean network size (by edge #):", mean(network.edge.sz), "\n")
cat(sum(network.edge.sz==1),"networks are single edge,", sum(network.edge.sz > 1), "are multi-edge.\n")
network.gene.sz <- unlist(lapply(1:max(gene.pair.table$NETWORK), function(net){length(unique(c(gene.pair.table[gene.pair.table$NETWORK==net,c("REG.GENE","TARGET.GENE")],recursive=T)))}))
cat("Range of network sizes (by gene #):", range(network.gene.sz), "\n")
cat("Mean network size (by gene #):", mean(network.gene.sz), "\n")

# Output the table
output.file <- paste0(my.args$OUTPUT, ".cis.trans.network.txt")
cat("Writing data on", nrow(gene.pair.table), "gene->gene network edges to:", output.file, "\n")
if(file.exists(output.file)) {
  cat("WARNING:", output.file, "exists already - overwriting now!\n")
}
write.table(gene.pair.table, output.file, sep="\t", row.names=F, quote=F)
cat("\n")

# Make a table summarizing the networks (number of genes, number of edges, density)
network.table <- data.frame(NETWORK=1:max(gene.pair.table$NETWORK), EDGE.SZ=as.integer(network.edge.sz), GENE.SZ=network.gene.sz)
network.table[,"DENSITY"] <- network.table$EDGE.SZ / (network.table$GENE.SZ * (network.table$GENE.SZ-1))
output.file <- paste0(my.args$OUTPUT, ".cis.trans.subnetworks.txt")
cat("Writing data on", nrow(network.table), "subnetworks to:", output.file, "\n")
if(file.exists(output.file)) {
  cat("WARNING:", output.file, "exists already - overwriting now!\n")
}
write.table(network.table, output.file, sep="\t", row.names=F, quote=F)

# Compute the number of targets and upstream regulators per gene, to determine likely "hub" genes
gene.table <- data.frame(GENE=unique(c(gene.pair.table$REG.GENE, gene.pair.table$TARGET.GENE)))
# Testing:
# gene.table <- data.frame(GENE=unique(c(gene.pair.table$REG.GENE, gene.pair.table$TARGET.GENE))[1:10])
gene.table$NETWORK <- foreach(gene=gene.table$GENE, .combine='c') %dopar% {
  unique(gene.pair.table[(gene.pair.table$REG.GENE==gene)|(gene.pair.table$TARGET.GENE==gene),"NETWORK"])
}
gene.table$DIRECT.TARGETS <- foreach(gene=gene.table$GENE, .combine='c') %dopar% {
  length(unique(gene.pair.table[gene.pair.table$REG.GENE==gene,"TARGET.GENE"]))
}
gene.table$TOTAL.TARGETS <- foreach(gene=gene.table$GENE, .combine='c') %dopar% {
  downstream.sz <- 0
  downstream <- c()
  repeat {
    downstream <- unique(gene.pair.table[gene.pair.table$REG.GENE %in% c(gene,downstream),"TARGET.GENE"])
    if(length(downstream) == downstream.sz) {
      break
    } else {
      downstream.sz <- length(downstream)
    }
  }
  return(length(setdiff(downstream, gene)))
}
gene.table$DIRECT.REGS <- foreach(gene=gene.table$GENE, .combine='c') %dopar% {
  length(unique(gene.pair.table[gene.pair.table$TARGET.GENE==gene,"REG.GENE"]))
}
gene.table$TOTAL.REGS <- foreach(gene=gene.table$GENE, .combine='c') %dopar% {
  upstream.sz <- 0
  upstream <- c()
  repeat {
    upstream <- unique(gene.pair.table[gene.pair.table$TARGET.GENE %in% c(gene,upstream),"REG.GENE"])
    if(length(upstream) == upstream.sz) {
      break
    } else {
      upstream.sz <- length(upstream)
    }
  }
  return(length(setdiff(upstream, gene)))
}

# Output this table
output.file <- paste0(my.args$OUTPUT, ".cis.trans.genes.txt")
cat("Writing data on", nrow(gene.table), "gene nodes to:", output.file, "\n")
if(file.exists(output.file)) {
  cat("WARNING:", output.file, "exists already - overwriting now!\n")
}
write.table(gene.table, output.file, sep="\t", row.names=F, quote=F)

cat("\nScript completed successfully!\n\n")
print(proc.time())
