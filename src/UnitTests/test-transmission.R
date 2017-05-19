p = list(
  num.acts = matrix(data=c(
                      100, 100, 100, 100,
                      100, 100, 100, 100,
                      100, 100,  10,  10,
                      100, 100,  10,   2), nrow=4, byrow=TRUE),
  p.condom = matrix(data=c(
                      0.02, 0.02, 0.02, 0.02,
                      0.02, 0.02, 0.02, 0.02,
                      0.02, 0.02, 0.15, 0.15,
                      0.02, 0.02, 0.15, 0.40), nrow=4, byrow=TRUE),
  breakup = 1.0,
  eff.mmc = 1.0 - 0.6,
  eff.condom = 1.0 - 0.9,
  infect.mult = c(13, 13, 2, 2, 4),
  infect.base = 0.002)

names = list(
  sex = c('M', 'W')
  )

ibm = read.table('output/ibm-test-transmission.txt', header=TRUE, comment.char='%')

for (g in 0:1) { ## donor sex
  for (h in 1:5) { ## stage
    for (kd in 0:3) { ## donor sexual activity level
      for (kr in 0:3) { ## recipient sexual activity level
        for (m in 0:1) { ## recipient MMC status
          ncdm = p$num.acts[kd+1,kr+1] * p$p.condom[kd+1,kr+1]
          nnon = p$num.acts[kd+1,kr+1] * (1.0 - p$p.condom[kd+1,kr+1])
          beta = p$infect.base * p$infect.mult[h]
          if (m == 0) {
            prob = 1 - ((1 - p$eff.condom * beta)^ncdm) * ((1 - beta)^nnon)
          } else {
            prob = 1 - ((1 - p$eff.condom * beta * p$eff.mmc)^ncdm) * ((1 - beta * p$eff.mmc)^nnon)
          }
          rate = p$breakup * prob / (1 - prob)

          indices = which((ibm$sp.sex == g) & (ibm$sp.level == kd) & (ibm$sp.stage == h) & (ibm$sn.level == kr) & (ibm$sn.mmc == m))
          if (length(indices) > 0) {
            ks.pval = ks.test(ibm$time[indices], "pexp", rate)$p

            xmin = min(ibm$time[indices])
            xmax = max(ibm$time[indices])
            rmin = 0
            rmax = 2^(ceiling(log(xmax, 2)))
            x = seq(rmin, rmax, (rmax - rmin) / 256)

            plot(c(rmin, rmax), c(0, 1), cex=0, axes=FALSE, ann=FALSE)
            lines(ecdf(ibm$time[indices]), lwd=2)
            lines(x, pexp(x, rate), col='#ff0000', lwd=2)
            axis(1)
            axis(2)
            title(main=sprintf("h=%d donor.g=%s, donor.k=%d recip.k=%d recip.mmc=%d\n(n=%d, p=%0.3f)\n", h, names$sex[g+1], kd, kr, m, length(indices), ks.pval))
            title(xlab="Transmission time")
            title(ylab="Cumulative probability")

            legend('bottomright', inset=0.025, box.lty=1, bg='#ffffff',
                   legend=c('Observed distribution (IBM)', 'Expected distribution'),
                   col=c('#000000', '#ff0000'), lty=1, lwd=2)

            readline("Press a key to see the next plot")
          } else {
            if (!(g == 0 | m == 1)) {
              ## We should not have observations for the case of an
              ## infected man (g==0) whose partner is circumcised,
              ## since we are only tracking heterosexual partnerships
              ## and MMC is for men only.
              cat(sprintf("No observations for %d %s %d %d %d\n", h, names$sex[g+1], kd, kr, m))
            }
          }
        }
      }
    }
  }
}
