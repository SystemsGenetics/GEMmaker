# Read in and log transform the matrix
ematrix = read.table('GEM.txt', header=TRUE, sep="\t", quote="")
samples = names(ematrix);
ematrix = log2(as.matrix(ematrix))

# Plot the samples
png('GEM_density_before.png')
colors <- rainbow(ncol(ematrix))
plot(density(ematrix[,1], na.rm=TRUE), xlab="log count")
for (i in 2:ncol(ematrix)) {
  if (length(which(is.na(ematrix[,i]))) == nrow(ematrix)) {
    next;
  }
  lines(density(ematrix[,i], na.rm=TRUE), col=colors[i])
}
dev.off()

# Plot the boxplots
png('GEM_boxplots_before.png')
boxplot(ematrix, las=2)
dev.off()


# Create a combined array of all values, use this for a Kolmogorov-Smirnov Outlier Test
g = ematrix[, 1]
for (i in 2:ncol(ematrix)) {
  g = c(g, ematrix[, i])
}

# Perform the KS test for each sample in the dataset:
ks_test = numeric()
for (i in 1:ncol(ematrix)) {
  if (length(which(is.na(ematrix[,i]))) == nrow(ematrix)) {
    ks_test[i] = 1;
    next;
  }
  ks_test[i] = ks.test(ematrix[, i], g)
}

# Convert results into a data frame:
ksdf = data.frame(samples, unlist(ks_test))
colnames(ksdf) = c('sample', 'ks_pvalue')
row.names(ksdf) = c(1:ncol(ematrix))

# Now that we have identified outliers we can use the following R code to remove them:
ks_th = 0.15
outliers = colnames(ematrix)[which(ksdf$ks_pvalue > ks_th)]
ematno = ematrix[,!(samples %in% outliers)]

# Plot again those that remain to see if the distribution functions are comparable
png('GEM_density_after.png')
colors <- rainbow(ncol(ematno))
plot(density(ematno[,1], na.rm=TRUE), xlab="log count")
for (i in 2:ncol(ematno)) {
  lines(density(ematno[,i], na.rm=TRUE), col=colors[i])
}
dev.off()

# Box Plots of just the non-outlier samples
png('GEM_boxplots_after.png')
boxplot(ematno, las=2)
dev.off()

# The following R code will output the expression matrix file
write.table(ksdf, file="GEM-ks_test_results.txt", append=FALSE, quote=FALSE, sep="\t", col.names=T)
write.table(ematno, file="GEM-log-no.txt", append=FALSE, quote=FALSE, sep="\t")
