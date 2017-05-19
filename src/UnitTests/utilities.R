ibm.risk.w = function(sim) {
  nwk0 = sim$nwk0b0 + sim$nwk0b1 + sim$nwk0b2 + sim$nwk0b3 + sim$nwk0b4 + sim$nwk0b5 + sim$nwk0b6 + sim$nwk0b7
  nwk1 = sim$nwk1b0 + sim$nwk1b1 + sim$nwk1b2 + sim$nwk1b3 + sim$nwk1b4 + sim$nwk1b5 + sim$nwk1b6 + sim$nwk1b7
  nwk2 = sim$nwk2b0 + sim$nwk2b1 + sim$nwk2b2 + sim$nwk2b3 + sim$nwk2b4 + sim$nwk2b5 + sim$nwk2b6 + sim$nwk2b7
  nwk3 = sim$nwk3b0 + sim$nwk3b1 + sim$nwk3b2 + sim$nwk3b3 + sim$nwk3b4 + sim$nwk3b5 + sim$nwk3b6 + sim$nwk3b7
  return(cbind(nwk0, nwk1, nwk2, nwk3))
}

ibm.risk.m = function(sim) {
  nmk0 = sim$nmk0b0 + sim$nmk0b1 + sim$nmk0b2 + sim$nmk0b3 + sim$nmk0b4 + sim$nmk0b5 + sim$nmk0b6 + sim$nmk0b7
  nmk1 = sim$nmk1b0 + sim$nmk1b1 + sim$nmk1b2 + sim$nmk1b3 + sim$nmk1b4 + sim$nmk1b5 + sim$nmk1b6 + sim$nmk1b7
  nmk2 = sim$nmk2b0 + sim$nmk2b1 + sim$nmk2b2 + sim$nmk2b3 + sim$nmk2b4 + sim$nmk2b5 + sim$nmk2b6 + sim$nmk2b7
  nmk3 = sim$nmk3b0 + sim$nmk3b1 + sim$nmk3b2 + sim$nmk3b3 + sim$nmk3b4 + sim$nmk3b5 + sim$nmk3b6 + sim$nmk3b7
  return(cbind(nmk0, nmk1, nmk2, nmk3))
}

ibm.prop.csw = function(sim) {
  nwk0 = sim$nwk0b0 + sim$nwk0b1 + sim$nwk0b2 + sim$nwk0b3 + sim$nwk0b4 + sim$nwk0b5 + sim$nwk0b6 + sim$nwk0b7
  nwk1 = sim$nwk1b0 + sim$nwk1b1 + sim$nwk1b2 + sim$nwk1b3 + sim$nwk1b4 + sim$nwk1b5 + sim$nwk1b6 + sim$nwk1b7
  nwk2 = sim$nwk2b0 + sim$nwk2b1 + sim$nwk2b2 + sim$nwk2b3 + sim$nwk2b4 + sim$nwk2b5 + sim$nwk2b6 + sim$nwk2b7
  nwk3 = sim$nwk3b0 + sim$nwk3b1 + sim$nwk3b2 + sim$nwk3b3 + sim$nwk3b4 + sim$nwk3b5 + sim$nwk3b6 + sim$nwk3b7
  return(nwk3 / (nwk0 + nwk1 + nwk2 + nwk3))
}

values = function(x) {
  mu = mean(x)
  vr = var(x)
  nm = length(x)
  ci = 1.96 * sqrt(vr / nm)
  return(c(mu, mu - ci, mu + ci))
}
