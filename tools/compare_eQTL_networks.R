#!/home/ljeveret/Tools/R-3.1.1/bin/Rscript
#
# LJE - 8/24/17
#
# compare_eQTL_networks.R
#
# Simple script to compare 2 eQTL networks in the format output by eQTL_network.R
#
# Usage:
# compare_eQTL_networks.R networkA.txt networkB.txt
#
# Both network files should be in the format:
# GENEX GENEY [Details]
# Meaning that GeneX is inferred to regulate GeneY, with details showing the SNPs involved, etc.
#

# --- INITIALIZATION --- #

# setwd("~/Projects/DGRP_Baseline_RNAseq_Post/")

options(stringsAsFactors=F)

usageStr="USAGE:\ncompare_eQTL_networks.R networkA.txt networkB.txt"
my.args <- commandArgs(trailingOnly=T)
# TEMP TESTING:
# my.args <- c("known_all_novel_genes/plink/ExprF.0.05.fdr.cis.trans.network.txt", "known_all_novel_genes/plink/results.OLD/ExprF.0.05.fdr.cis.trans.network.txt")

if(length(my.args) != 2) {
  stop("Requires exactly two arguments. ", usageStr)
}

# Parse param names and values
net.A.file <- my.args[1]
net.B.file <- my.args[2]

if(!file.exists(net.A.file)) {
  stop(net.A.file, " does not exist.")
}

if(!file.exists(net.B.file)) {
  stop(net.B.file, " does not exist.")
}

# Load both files
cat("Loading Network A =", net.A.file, "\n")
net.A <- read.table(net.A.file, header=T, sep="\t")
cat("Loading Network B =", net.B.file, "\n")
net.B <- read.table(net.B.file, header=T, sep="\t")

cat("Loaded", nrow(net.A), "edges for Network A\n")
cat("Loaded", nrow(net.B), "edges for Network B\n")
cat("\n")

# Extract list of nodes (reg, target, and union of both)
net.A.reg <- unique(net.A[,1])
net.A.target <- unique(net.A[,2])
net.A.nodes <- union(net.A.reg, net.A.target)

net.B.reg <- unique(net.B[,1])
net.B.target <- unique(net.B[,2])
net.B.nodes <- union(net.B.reg, net.B.target)

shared.reg <- intersect(net.A.reg, net.B.reg)
shared.target <- intersect(net.A.target, net.B.target)
shared.nodes <- intersect(net.A.nodes, net.B.nodes)

all.reg <- union(net.A.reg, net.B.reg)
all.target <- union(net.A.target, net.B.target)
all.nodes <- union(net.A.nodes, net.B.nodes)

net.A.reg.only <- setdiff(net.A.reg, net.B.reg)
net.A.target.only <- setdiff(net.A.target, net.B.target)
net.A.node.only <- setdiff(net.A.nodes, net.B.nodes)

net.B.reg.only <- setdiff(net.B.reg, net.A.reg)
net.B.target.only <- setdiff(net.B.target, net.A.target)
net.B.node.only <- setdiff(net.B.nodes, net.A.nodes)

cat("Network A contains", length(net.A.nodes), "nodes.\n")
cat("Network B contains", length(net.B.nodes), "nodes.\n")
cat(length(shared.nodes), " nodes (", round(length(shared.nodes)*100/length(all.nodes)), "%) are in both networks.\n", sep="")
cat(length(net.A.node.only), " nodes (", round(length(net.A.node.only)*100/length(all.nodes)), "%) are in Network A ONLY.\n", sep="")
cat(length(net.B.node.only), " nodes (", round(length(net.B.node.only)*100/length(all.nodes)), "%) are in Network B ONLY.\n", sep="")
cat("\n")

cat("Network A contains", length(net.A.reg), "regulators.\n")
cat("Network B contains", length(net.B.reg), "regulators.\n")
cat(length(shared.reg), " regulators (", round(length(shared.reg)*100/length(all.reg)), "%) are in both networks.\n", sep="")
cat(length(net.A.reg.only), " regulators (", round(length(net.A.reg.only)*100/length(all.reg)), "%) are in Network A ONLY.\n", sep="")
cat(length(net.B.reg.only), " regulators (", round(length(net.B.reg.only)*100/length(all.reg)), "%) are in Network B ONLY.\n", sep="")
cat("\n")

cat("Network A contains", length(net.A.target), "targets.\n")
cat("Network B contains", length(net.B.target), "targets.\n")
cat(length(shared.target), " targets (", round(length(shared.target)*100/length(all.target)), "%) are in both networks.\n", sep="")
cat(length(net.A.target.only), " targets (", round(length(net.A.target.only)*100/length(all.target)), "%) are in Network A ONLY.\n", sep="")
cat(length(net.B.target.only), " targets (", round(length(net.B.target.only)*100/length(all.target)), "%) are in Network B ONLY.\n", sep="")
cat("\n")

net.A.edges <- unique(paste(net.A[,1], net.A[,2], sep=":"))
net.B.edges <- unique(paste(net.B[,1], net.B[,2], sep=":"))
stopifnot(length(net.A.edges)==nrow(net.A))
stopifnot(length(net.B.edges)==nrow(net.B))

shared.edges <- intersect(net.A.edges, net.B.edges)
all.edges <- union(net.A.edges, net.B.edges)
net.A.edges.only <- setdiff(net.A.edges, net.B.edges)
net.B.edges.only <- setdiff(net.B.edges, net.A.edges)

cat("Network A contains", length(net.A.edges), "edges.\n")
cat("Network B contains", length(net.B.edges), "edges.\n")
cat(length(shared.edges), " edges (", round(length(shared.edges)*100/length(all.edges)), "%) are in both networks.\n", sep="")
cat(length(net.A.edges.only), " edges (", round(length(net.A.edges.only)*100/length(all.edges)), "%) are in network A ONLY.\n", sep="")
cat(length(net.B.edges.only), " edges (", round(length(net.B.edges.only)*100/length(all.edges)), "%) are in network B ONLY.\n", sep="")

cat("\nScript completed successfully!\n")
