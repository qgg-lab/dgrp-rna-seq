# ========================
# = analyze blat results =
# ========================

blat <- read.table("all.blast8", header = FALSE, as.is = TRUE)
blat.top2 <- t(sapply(split(blat$V12, blat$V1), function(x) { return(sort(c(0, x), decreasing = T)[1:2]) }))

save(blat.top2, file = "blat.top2.RData")

ntr.id <- read.table("ntr.id", header = F, as.is = T)
bad.id <- read.table("../bad.gene.txt", header = F, as.is = T)[, 1]
rownames(ntr.id) <- ntr.id[, 1]

blat.top2 <- cbind(blat.top2, ifelse(grepl("TCON", rownames(blat.top2)), 1, 0))
blat.top2 <- cbind(blat.top2, ifelse((ntr.id[rownames(blat.top2), 2]) %in% bad.id, 1, 0))
blat.top2 <- cbind(blat.top2, blat.top2[, 1]-blat.top2[, 2])

# make graph

file.width = 60
cairo_pdf(file = "bitscore.pdf", width = file.width/25.4, height = file.width/25.4*0.8, family = "Arial")
par(las = 1, tcl = -0.2, mar = c(1.5, 2, 0.5, 0), ps = 7, lwd = 0.5, xpd = TRUE)

plot(c(0, 5), c(0, 1.4), type = "n", xlab = "", ylab = "", axes = FALSE)

# known gene
lines(density(log10(blat.top2[blat.top2[, 3] == 0, 5] + 1), bw = 0.05, from = -0.1, to = 5, n = 510), col = "black", lwd = 2)

# ntr filtered
lines(density(log10(blat.top2[blat.top2[, 3] == 1 & blat.top2[, 4] == 1, 5] + 1), bw = 0.05, from = -0.1, to = 5, n = 510), col = "red", lwd = 2)

# ntr retained
lines(density(log10(blat.top2[blat.top2[, 3] == 1 & blat.top2[, 4] == 0, 5] + 1), bw = 0.05, from = -0.1, to = 5, n = 510), col = "blue", lwd = 2)


axis(side = 1, lwd = 0.5, cex.axis = 6/par("ps")/par("cex"), mgp = c(2.5, 0, 0))
axis(side = 2, lwd = 0.5, cex.axis = 6/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
box(bty = "l", lwd = 0.5)
title(xlab = expression(paste("log"[10], "(", Delta, "Bitscore + 1)")), mgp = c(0.65, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Density", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("topleft", lwd = 2, bty = "n", legend = c("Known genes", "NTR filtered", "NTR retained"), col = c("black", "red", "blue"), x.intersp = 0.5, y.intersp = 0.5)

dev.off()

