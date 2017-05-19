path='output'

ode = read.table(sprintf('%s/ode-test-partnership.txt', path), header=TRUE, comment.char='%')
ibm = list()
for (k in 1:100) {
  ibm[[k]] = read.table(sprintf('%s/ibm-test-partnership-%02d.txt', path, k-1), header=TRUE, comment.char='%')
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
pdf("trend-popsize.pdf", h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,8e3), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ode$year, ode$n, col='#000000', lwd=1.5, lty=1)
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#ff0000', lty=1, lwd=1.5)
axis(1)
axis(2)
title(xlab="Year")
title(ylab="Overall population size")
dev.off()

ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X$couples)})), ncol=length(ibm))
pdf("trend-couples.pdf", h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,150), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ode$year, ode$couples, col='#000000', lwd=1.5, lty=1)
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#ff0000', lty=1, lwd=1.5)
axis(1)
axis(2)
title(xlab="Year")
title(ylab="Number of couples")
legend('topleft', inset=0.025, box.lty=1, bg='#ffffff',
       legend=c('ODE', 'IBM (average)', 'IBM (single runs)'),
       col=c('#000000', '#ff0000', '#666666'),
       lty=1, lwd=c(1.5, 1.5, 1))
dev.off()

print(c(ode$couples[nrow(ode)], values(ibm.X[nrow(ibm.X),])))

## different age, different risk
ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X$nw21m42)})), ncol=length(ibm))
pdf("trend-nw21m42.pdf", h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,6), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ode$year, ode$nw21m42, col='#000000', lwd=1.5, lty=1)
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#ff0000', lty=1, lwd=1.5)
axis(1)
axis(2)
title(xlab="Year")
title(ylab="Number of couples")
dev.off()

## same age, different risk 
ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X$nw21m22)})), ncol=length(ibm))
pdf("trend-nw21m22.pdf", h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,25), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ode$year, ode$nw21m22, col='#000000', lwd=1.5, lty=1)
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#ff0000', lty=1, lwd=1.5)
axis(1)
axis(2)
title(xlab="Year")
title(ylab="Number of couples")
dev.off()

## disparate age, different risk 
ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X$nw21m32)})), ncol=length(ibm))
pdf("trend-nw21m32.pdf", h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,20), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ode$year, ode$nw21m32, col='#000000', lwd=1.5, lty=1)
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#ff0000', lty=1, lwd=1.5)
axis(1)
axis(2)
title(xlab="Year")
title(ylab="Number of couples")
dev.off()

## different age, same risk 
ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X$nw21m41)})), ncol=length(ibm))
pdf("trend-nw21m41.pdf", h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,10), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ode$year, ode$nw21m41, col='#000000', lwd=1.5, lty=1)
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#ff0000', lty=1, lwd=1.5)
axis(1)
axis(2)
title(xlab="Year")
title(ylab="Number of couples")
dev.off()

## same age, same risk 
ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X$nw21m21)})), ncol=length(ibm))
pdf("trend-nw21m21.pdf", h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,30), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ode$year, ode$nw21m21, col='#000000', lwd=1.5, lty=1)
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#ff0000', lty=1, lwd=1.5)
axis(1)
axis(2)
title(xlab="Year")
title(ylab="Number of couples")
dev.off()

## disparate age, same risk 
ibm.X = matrix(unlist(lapply(ibm, function(X) {return(X$nw21m31)})), ncol=length(ibm))
pdf("trend-nw21m31.pdf", h=5)
par(las=1, mar=c(5,5,4,2)+0.1)
plot(range(ode$year), c(0,30), cex=0, ann=FALSE, axes=FALSE)
for (k in 1:ncol(ibm.X)) {lines(ibm[[k]]$year + 1978, ibm.X[,k], col=cset[k])}
lines(ode$year, ode$nw21m31, col='#000000', lwd=1.5, lty=1)
lines(ibm[[1]]$year + 1978, apply(ibm.X, 1, mean), col='#ff0000', lty=1, lwd=1.5)
axis(1)
axis(2)
title(xlab="Year")
title(ylab="Number of couples")
dev.off()
