ode = read.table('output/ode-test-demographer.txt', header=TRUE, comment.char='%')
ibm = list()
for (k in 1:100) {
  ibm[[k]] = read.table(sprintf('output/ibm-test-demographer-%02d.txt', k-1), header=TRUE, comment.char='%')
}

values = function(x) {
  mu = mean(x)
  vr = var(x)
  nm = length(x)
  ci = 1.96 * sqrt(vr / nm)
  return(c(mu, mu - ci, mu + ci))
}

cset=grey(runif(length(ibm), 0.6, 0.8))

ibm.X = matrix(unlist(lapply(ibm, function(X) {return(rowSums(X[,2:7]))})), ncol=length(ibm))
pdf('Figures/trend-popsize.pdf', h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,8e3), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#000000', lty=1, lwd=2)
lines(ode$year, rowSums(ode[,2:8]), col='#ff0000', lwd=2, lty=1)
axis(1)
axis(2)
title(xlab="Year")
title(ylab="Overall population size")
legend('topleft', inset=0.025, box.lty=1, bg='#ffffff',
       legend=c('ODE', 'IBM (average)', 'IBM (single run)'),
       col=c('#ff0000', '#000000', cset[1]),
       lwd=c(2,2,1))
dev.off()

bands=c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54")
risks=c("Least", "Low", "Medium", "High")

for (ki in 0:3) {
  for (ai in 0:7) {
    cname = sprintf("nwk%db%d", ki, ai)
    i.ibm = which(colnames(ibm[[1]]) == cname)
    i.ode = which(colnames(ode) == cname)
    ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X[,i.ibm])})), ncol=length(ibm))
    pdf(sprintf("Figures/trend-%s.pdf", cname), h=5)
    par(las=1, mar=c(5,5,4,2)+0.1)
    plot(range(ode$year), c(0,max(ibm.X)), cex=0, ann=FALSE, axes=FALSE)
    for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
    lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#ff0000', lty=1, lwd=2)
    lines(ode$year, ode[,i.ode], col='#000000', lwd=2, lty=1)
    axis(1)
    axis(2)
    title(ylab=sprintf("%s-risk women aged %s", risks[ki+1], bands[ai+1]))
    title(xlab="Year")
    legend('topleft', inset=0.025, box.lty=1, bg='#ffffff',
           legend=c('ODE', 'IBM (average)', 'IBM (single run)'),
           col=c('#ff0000', '#000000', cset[1]),
           lwd=c(2,2,1))
    dev.off()
  }
}

for (ki in 0:3) {
  for (ai in 0:7) {
    cname = sprintf("nmk%db%d", ki, ai)
    i.ibm = which(colnames(ibm[[1]]) == cname)
    i.ode = which(colnames(ode) == cname)
    ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X[,i.ibm])})), ncol=length(ibm))
    pdf(sprintf("Figures/trend-%s.pdf", cname), h=5)
    par(las=1, mar=c(5,5,4,2)+0.1)
    plot(range(ode$year), c(0,max(ibm.X)), cex=0, ann=FALSE, axes=FALSE)
    for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
    lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#ff0000', lty=1, lwd=2)
    lines(ode$year, ode[,i.ode], col='#000000', lwd=2, lty=1)
    axis(1)
    axis(2)
    title(ylab=sprintf("%s-risk men aged %s", risks[ki+1], bands[ai+1]))
    title(xlab="Year")
    legend('topleft', inset=0.025, box.lty=1, bg='#ffffff',
           legend=c('ODE', 'IBM (average)', 'IBM (single run)'),
           col=c('#ff0000', '#000000', cset[1]),
           lwd=c(2,2,1))
    dev.off()
  }
}
