source('utilities.R')

ode = read.table('output/ode-test-progression.txt', header=TRUE, comment.char='%')
ibm = list()
for (k in 1:100) {
  ibm[[k]] = read.table(sprintf('output/ibm-test-progression-%02d.txt', k-1), header=TRUE, comment.char='%')
}

cset=grey(runif(length(ibm), 0.6, 0.8))

ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X$P1)})), ncol=length(ibm))
pdf("Figures/trend-P1.pdf", h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,max(ibm.X)), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#000000', lty=1, lwd=1.5)
lines(ode$year, ode$P1, col='#ff0000', lwd=1.5, lty=1)
axis(1)
axis(2)
title(ylab="Acute window")
title(xlab="Year")
legend('topright', inset=0.025, box.lty=1, bg='#ffffff',
       legend=c('ODE', 'IBM (average)', 'IBM (single run)'),
       col=c('#ff0000', '#000000', grey(0.7)), lty=1, lwd=c(2,2,1))
dev.off()

ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X$P2)})), ncol=length(ibm))
pdf("Figures/trend-P2.pdf", h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,max(ibm.X)), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#000000', lty=1, lwd=1.5)
lines(ode$year, ode$P2, col='#ff0000', lwd=1.5, lty=1)
axis(1)
axis(2)
title(ylab="Acute established")
title(xlab="Year")
legend('topright', inset=0.025, box.lty=1, bg='#ffffff',
       legend=c('ODE', 'IBM (average)', 'IBM (single run)'),
       col=c('#ff0000', '#000000', grey(0.7)), lty=1, lwd=c(2,2,1))
dev.off()

ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X$L1)})), ncol=length(ibm))
pdf("Figures/trend-L1.pdf", h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,max(ibm.X)), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#000000', lty=1, lwd=1.5)
lines(ode$year, ode$L2, col='#ff0000', lwd=1.5, lty=1)
axis(1)
axis(2)
title(ylab="Early chronic")
title(xlab="Year")
legend('topright', inset=0.025, box.lty=1, bg='#ffffff',
       legend=c('ODE', 'IBM (average)', 'IBM (single run)'),
       col=c('#ff0000', '#000000', grey(0.7)), lty=1, lwd=c(2,2,1))
dev.off()

ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X$L2)})), ncol=length(ibm))
pdf("Figures/trend-L2.pdf", h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,max(ibm.X)), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#000000', lty=1, lwd=1.5)
lines(ode$year, ode$L3, col='#ff0000', lwd=1.5, lty=1)
axis(1)
axis(2)
title(ylab="Late chronic")
title(xlab="Year")
legend('topright', inset=0.025, box.lty=1, bg='#ffffff',
       legend=c('ODE', 'IBM (average)', 'IBM (single run)'),
       col=c('#ff0000', '#000000', grey(0.7)), lty=1, lwd=c(2,2,1))
dev.off()

ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X$A1)})), ncol=length(ibm))
pdf("Figures/trend-A1.pdf", h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,max(ibm.X)), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#000000', lty=1, lwd=1.5)
lines(ode$year, ode$A1, col='#ff0000', lwd=1.5, lty=1)
axis(1)
axis(2)
title(ylab="AIDS")
title(xlab="Year")
legend('topright', inset=0.025, box.lty=1, bg='#ffffff',
       legend=c('ODE', 'IBM (average)', 'IBM (single run)'),
       col=c('#ff0000', '#000000', grey(0.7)), lty=1, lwd=c(2,2,1))
dev.off()
