# =====================
# = summarize results =
# =====================

# read experiment info
# ============================================================

exp <- read.csv("tabula-TableS2.csv", header = T, as.is = T)

# try a sample
# ============================================================

i = 1
this.exp.id <- exp[i, 4]
this.quant <- read.table(paste("quant/", this.exp.id, ".quant/abundance.tsv", sep =""), header = TRUE, as.is = T)
tpm.matrix <- matrix(NA, ncol = nrow(exp), nrow = nrow(this.quant))
gene <- this.quant[, 1]

for (i in 1:nrow(exp)) {
	
	this.exp.id <- exp[i, 4]
	this.quant <- read.table(paste("quant/", this.exp.id, ".quant/abundance.tsv", sep =""), header = TRUE, as.is = T)
	tpm.matrix[, i] <- this.quant[, 5]
	
}

save(exp, tpm.matrix, gene, file = "dgrp-cell-line-exp.RData")

# IDs
ntr.id <- read.table("ntr.id", header = F, as.is = T)
bad.id <- read.table("../bad.gene.txt", header = F, as.is = T)[, 1]
rownames(ntr.id) <- ntr.id[, 1]

# ntr versus known
ntr.idx <- ifelse(grepl("TCON", gene), 1, 0)
bad.idx <- ifelse(ntr.id[gene, 2] %in% bad.id, 1, 0)

# find maximum expression
tpm.max <- apply(tpm.matrix, 1, max)
tpm.median <- apply(tpm.matrix, 1, median)

# make plots
# ============================================================

file.width = 100
cairo_pdf(file = "tpm.pdf", width = file.width/25.4, height = file.width/25.4*0.55, family = "Arial")
par(las = 1, tcl = -0.2, mar = c(1.5, 2, 1.5, 0), ps = 7, lwd = 0.5, xpd = TRUE, mfrow = c(1, 2))


plot(c(-4.1, 5), c(0, 3), type = "n", xlab = "", ylab = "", axes = FALSE)

# known gene
lines(density(log10(tpm.median + 0.0001), bw = 0.10, from = -4.2, to = 5, n = 512), col = "black", lwd = 2)

# ntr filtered
lines(density(log10(tpm.median[ntr.idx == 1 & bad.idx == 1] + 0.0001), bw = 0.10, from = -4.2, to = 5, n = 512), col = "red", lwd = 2)

# ntr retained
lines(density(log10(tpm.median[ntr.idx == 1 & bad.idx == 0] + 0.0001), bw = 0.10, from = -4.2, to = 5, n = 512), col = "blue", lwd = 2)


axis(side = 1, lwd = 0.5, cex.axis = 6/par("ps")/par("cex"), mgp = c(2.5, 0, 0))
axis(side = 2, lwd = 0.5, cex.axis = 6/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
box(bty = "l", lwd = 0.5)
title(xlab = expression(paste("log"[10], "(TPM + 0.0001)")), mgp = c(0.65, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Density", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "Median among cell lines")

legend("topleft", lwd = 2, bty = "n", legend = c("Known genes", "NTR filtered", "NTR retained"), col = c("black", "red", "blue"), x.intersp = 0.5, y.intersp = 0.5)

plot(c(-4.1, 5), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)

# known gene
lines(density(log10(tpm.max + 0.0001), bw = 0.10, from = -4.2, to = 5, n = 512), col = "black", lwd = 2)

# ntr filtered
lines(density(log10(tpm.max[ntr.idx == 1 & bad.idx == 1] + 0.0001), bw = 0.10, from = -4.2, to = 5, n = 512), col = "red", lwd = 2)

# ntr retained
lines(density(log10(tpm.max[ntr.idx == 1 & bad.idx == 0] + 0.0001), bw = 0.10, from = -4.2, to = 5, n = 512), col = "blue", lwd = 2)


axis(side = 1, lwd = 0.5, cex.axis = 6/par("ps")/par("cex"), mgp = c(2.5, 0, 0))
axis(side = 2, lwd = 0.5, cex.axis = 6/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
box(bty = "l", lwd = 0.5)
title(xlab = expression(paste("log"[10], "(TPM + 0.0001)")), mgp = c(0.65, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Density", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "Max among cell lines")
legend("topleft", lwd = 2, bty = "n", legend = c("Known genes", "NTR filtered", "NTR retained"), col = c("black", "red", "blue"), x.intersp = 0.5, y.intersp = 0.5)



dev.off()



