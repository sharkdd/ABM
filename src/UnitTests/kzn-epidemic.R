## Script used to visualize output from kzn-epidemic.sh

ibm = list()
ibm.files = list.files("output", pattern="ibm-kzn-epi-[0-9]+.txt")
for (k in 1:length(ibm.files)) {
  ibm[[k]] = read.table(sprintf("output/%s", ibm.files[k]), header=TRUE, comment.char='%')
}
ode = read.table('output/ode-kzn-epi.txt', header=TRUE, comment.char='%')

visualize = function(ibm.year, ibm.data, ode.year, ode.data, descr) {
  ibm.avg = apply(ibm.data, 1, mean)
  cset=grey(runif(ncol(ibm.data), 0.6, 0.8))
  ymin = min(cbind(ode.data, ibm.data))
  ymax = max(cbind(ode.data, ibm.data))

  if (length(dev.list()) == 0) {quartz(h=5)}
  par(las=1, mar=c(5,5,4,2)+0.1)
  plot(range(ode.year), c(ymin, ymax), cex=0, ann=FALSE, axes=FALSE)
  for (k in 1:ncol(ibm.data)) {
    lines(ibm.year[,k], ibm.data[,k], col=cset[k])
  }
  lines(ibm.year[,1], ibm.avg, col='#000000', lwd=2)
  lines(ode.year, ode.data, col='#ff0000', lwd=2)
  axis(1)
  axis(2)
  legend('topleft', inset=0.025, box.lty=1, bg='#ffffff',
         legend=c('ODE', 'IBM (average)', 'IBM (single run)'),
         col=c('#ff0000', '#000000', grey(0.7)), lty=1, lwd=c(2,2,1))
  title(xlab="Year")
  title(ylab=descr)
}

## The IBM starts reports simulation time starting from zero, while
## the ODE reports simulation time starting from 1978.
ibm.year = data.frame(sapply(ibm, function(sim) {return(sim$year + 1978)}))

ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$X)}))
visualize(ibm.year, ibm.data, ode$year, ode$X, "Susceptible")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$P1)}))
visualize(ibm.year, ibm.data, ode$year, ode$P1, "Window acute infection")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$P2)}))
visualize(ibm.year, ibm.data, ode$year, ode$P2, "Established acute infection")

## The KZN IBM represents two stages of chronic infection (L1 and
## L2). The KZN ODE has support for three (L1, L2 and L3), but only
## uses the latter two.
readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$L1)}))
visualize(ibm.year, ibm.data, ode$year, ode$L2, "Chronic infection, CD4>350")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$L2)}))
visualize(ibm.year, ibm.data, ode$year, ode$L3, "Chronic infection, CD4 201-350")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$A1)}))
visualize(ibm.year, ibm.data, ode$year, ode$A1, "AIDS")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$ART.L2 + sim$ART.A1)}))
visualize(ibm.year, ibm.data, ode$year, ode$art.L3 + ode$art.A1, "Number on ART")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$DR)}))
visualize(ibm.year, ibm.data, ode$year, ode$DR, "Cases of overall drug resistance")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$mmc.X + sim$mmc.Y)}))
visualize(ibm.year, ibm.data, ode$year, ode$mmc.X + ode$mmc.Y, "Circumcised men")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$nwk0b0 + sim$nwk0b1 + sim$nwk0b2 + sim$nwk0b3 + sim$nwk0b4 + sim$nwk0b5 + sim$nwk0b6 + sim$nwk0b7)}))
visualize(ibm.year, ibm.data, ode$year, ode$nwk0b0 + ode$nwk0b1 + ode$nwk0b2 + ode$nwk0b3 + ode$nwk0b4 + ode$nwk0b5 + ode$nwk0b6 + ode$nwk0b7, "Women, least activity level")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$nwk1b0 + sim$nwk1b1 + sim$nwk1b2 + sim$nwk1b3 + sim$nwk1b4 + sim$nwk1b5 + sim$nwk1b6 + sim$nwk1b7)}))
visualize(ibm.year, ibm.data, ode$year, ode$nwk1b0 + ode$nwk1b1 + ode$nwk1b2 + ode$nwk1b3 + ode$nwk1b4 + ode$nwk1b5 + ode$nwk1b6 + ode$nwk1b7, "Women, low activity level")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$nwk2b0 + sim$nwk2b1 + sim$nwk2b2 + sim$nwk2b3 + sim$nwk2b4 + sim$nwk2b5 + sim$nwk2b6 + sim$nwk2b7)}))
visualize(ibm.year, ibm.data, ode$year, ode$nwk2b0 + ode$nwk2b1 + ode$nwk2b2 + ode$nwk2b3 + ode$nwk2b4 + ode$nwk2b5 + ode$nwk2b6 + ode$nwk2b7, "Women, medium activity level")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$nwk3b0 + sim$nwk3b1 + sim$nwk3b2 + sim$nwk3b3 + sim$nwk3b4 + sim$nwk3b5 + sim$nwk3b6 + sim$nwk3b7)}))
visualize(ibm.year, ibm.data, ode$year, ode$nwk3b0 + ode$nwk3b1 + ode$nwk3b2 + ode$nwk3b3 + ode$nwk3b4 + ode$nwk3b5 + ode$nwk3b6 + ode$nwk3b7, "Women, high activity level")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$nmk0b0 + sim$nmk0b1 + sim$nmk0b2 + sim$nmk0b3 + sim$nmk0b4 + sim$nmk0b5 + sim$nmk0b6 + sim$nmk0b7)}))
visualize(ibm.year, ibm.data, ode$year, ode$nmk0b0 + ode$nmk0b1 + ode$nmk0b2 + ode$nmk0b3 + ode$nmk0b4 + ode$nmk0b5 + ode$nmk0b6 + ode$nmk0b7, "Women, least activity level")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$nmk1b0 + sim$nmk1b1 + sim$nmk1b2 + sim$nmk1b3 + sim$nmk1b4 + sim$nmk1b5 + sim$nmk1b6 + sim$nmk1b7)}))
visualize(ibm.year, ibm.data, ode$year, ode$nmk1b0 + ode$nmk1b1 + ode$nmk1b2 + ode$nmk1b3 + ode$nmk1b4 + ode$nmk1b5 + ode$nmk1b6 + ode$nmk1b7, "Women, low activity level")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$nmk2b0 + sim$nmk2b1 + sim$nmk2b2 + sim$nmk2b3 + sim$nmk2b4 + sim$nmk2b5 + sim$nmk2b6 + sim$nmk2b7)}))
visualize(ibm.year, ibm.data, ode$year, ode$nmk2b0 + ode$nmk2b1 + ode$nmk2b2 + ode$nmk2b3 + ode$nmk2b4 + ode$nmk2b5 + ode$nmk2b6 + ode$nmk2b7, "Women, medium activity level")

readline("Press a key to see the next plot")
ibm.data = data.frame(sapply(ibm, function(sim) {return(sim$nmk3b0 + sim$nmk3b1 + sim$nmk3b2 + sim$nmk3b3 + sim$nmk3b4 + sim$nmk3b5 + sim$nmk3b6 + sim$nmk3b7)}))
visualize(ibm.year, ibm.data, ode$year, ode$nmk3b0 + ode$nmk3b1 + ode$nmk3b2 + ode$nmk3b3 + ode$nmk3b4 + ode$nmk3b5 + ode$nmk3b6 + ode$nmk3b7, "Women, high activity level")
