library(rethinking)
data(rugged)
d <- rugged

# make log version of outcome
d$log_gdp <- log(d$rgdppc_2000)

# extract countries with GDP data
dd <- d[complete.cases(d$rgdppc_2000), ]

# rescale variables
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)

r_bar <- mean(dd$rugged_std)

m8.1 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a + b * (rugged_std - r_bar),
    a ~ dnorm(1, 1),
    b ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=dd
)

set.seed(7)
prior <- extract.prior(m8.1)

# set up the plot dimensions
plot(NULL, xlim=c(0, 1), ylim=c(0.5, 1.5), xlab="ruggedness", ylab="log GDP")
abline(h=min(dd$log_gdp_std), lty=2)
abline(h=max(dd$log_gdp_std), lty=2)

# draw 50 lines from the prior
rugged_seq <- seq(from=-0.1, to=1.1, length.out=30)
mu <- link(m8.1, post=prior, data=data.frame(rugged_std=rugged_seq))
for(i in 1:50)
  lines(rugged_seq, mu[i, ], col=col.alpha("black", 0.3))


m8.1 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a + b * (rugged_std - r_bar),
    a ~ dnorm(1, 0.1),
    b ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ), data=dd
)

precis(m8.1)

# make variable to index Africa (1) or not Africa (2)
dd$cid <- ifelse(dd$cont_africa == 1, 1, 2)

m8.2 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b * (rugged_std - r_bar),
    a[cid] ~ dnorm(1, 0.1),
    b ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ), data=dd
)

compare(m8.1, m8.2)

rugged.seq <- seq(from=-0.1, to=1.1, length.out=30)
mu.notAfrica <- link(m8.2, data=data.frame(cid=2, rugged_std=rugged.seq))
mu.africa <- link(m8.2, data=data.frame(cid=1, rugged_std=rugged.seq))
mu.notAfrica.mu <- apply(mu.notAfrica, 2, mean)
mu.notAfrica.pi <- apply(mu.notAfrica, 2, PI, prob=0.97)
mu.africa.mu <- apply(mu.africa, 2, mean)
mu.africa.pi <- apply(mu.africa, 2, PI, prob=.97)

## TODO: Add Plotting code (8.4)

m8.3 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b[cid] * (rugged_std - r_bar),
    a[cid] ~ dnorm(1, 0.1),
    b[cid] ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ), data=dd
)

precis(m8.3, depth=2)

compare(m8.1, m8.2, m8.3, func=PSIS)
plot(PSIS(m8.3, pointwise=TRUE)$k)

# plot Africa - cid=1
d.A1 <- dd[dd$cid==1, ]
plot(d.A1$rugged_std, d.A1$log_gdp_std, pch=16, col=rangi2, xlab="ruggedness (standardized)",
     ylab="log GDP (as proportion of mean)", xlim=c(0, 1))
mu <- link(m8.3, data=data.frame(cid=1, rugged_std=rugged_seq))
mu.mean <- apply(mu, 2, mean)
mu.ci <- apply(mu, 2, PI, prob=.97)
lines(rugged_seq, mu.mean, lwd=2)
shade(mu.ci, rugged_seq, col=col.alpha(rangi2, 0.3))
mtext("African nations")


# plot non Africa - cid=2
d.A0 <- dd[dd$cid==2, ]
plot(d.A0$rugged_std, d.A0$log_gdp_std, pch=1, col="black", xlab="ruggedness (standardized)",
     ylab="log GDP (as proportion of mean)", xlim=c(0, 1))
mu <- link(m8.3, data=data.frame(cid=2, rugged_std=rugged_seq))
mu.mean <- apply(mu, 2, mean)
mu.ci <- apply(mu, 2, PI, prob=.97)
lines(rugged_seq, mu.mean, lwd=2)
shade(mu.ci, rugged_seq)
mtext("Non-African nations")


rugged_seq <- seq(from=-.2, to=1.2, length.out=30)
muA <- link(m8.3, data=data.frame(cid=1, rugged_std=rugged_seq))
muN <- link(m8.3, data=data.frame(cid=2, rugged_std=rugged_seq))
delta <- muA - muN

delta.mean <- apply(delta, 2, mean)
delta.ci <- apply(delta, 2, PI, prob=.97)


plot(NULL, xlab="ruggedness", ylab="Expected difference African and non African nation", xlim=range(rugged_seq), ylim=c(-0.3, 0.2))
lines(rugged_seq, delta.mean, lwd=2)
shade(delta.ci, rugged_seq)


rm(list=ls())

### 8.3

data(tulips)
d <- tulips
str(d)

d$blooms_std <- d$blooms / max(d$blooms)
d$water_cent <- d$water - mean(d$water)
d$shade_cent <- d$shade - mean(d$shade)

m8.4 <- quap(
  alist(
    blooms_std ~ dnorm(mu, sigma),
    mu <- a + bw * water_cent + bs * shade_cent,
    a ~ dnorm(.5, .25),
    bw ~ dnorm(0, .25),
    bs ~ dnorm(0, .25),
    sigma ~ dexp(1)
  ), data=d
)

m8.5 <- quap(
  alist(
    blooms_std ~ dnorm(mu, sigma),
    mu <- a + bw * water_cent + bs * shade_cent + bws * water_cent * shade_cent,
    a ~ dnorm(.5, .25),
    bw ~ dnorm(0, .25),
    bs ~ dnorm(0, .25),
    bws ~ dnorm(0, .25),
    sigma ~ dexp(1)
  ), data=d
)

# Triptych plot

par(mfrow=c(2, 3))
for(model in c(m8.4, m8.5)) {
  for(s in -1:1) {
    idx <- which(d$shade_cent == s)
    plot(d$water_cent[idx], d$blooms_std[idx], xlim=c(-1, 1), ylim=c(0, 1), 
         xlab="water", ylab="blooms", pch=16, col=rangi2)
    mu <- link(model, data=data.frame(shade_cent=s, water_cent=-1:1))
    for(i in 1:20) lines(-1:1,  mu[i, ], col=col.alpha("black", .3))
  }
}

par(mfrow=c(2, 3))
for(model in c(m8.4, m8.5)) {
  prior <- extract.prior(model)
  for(s in -1:1) {
    idx <- which(d$shade_cent == s)
    plot(d$water_cent[idx], d$blooms_std[idx], xlim=c(-1, 1), ylim=c(0, 1), 
         xlab="water", ylab="blooms", pch=16, col=rangi2)
    mu <- link(model, post=prior, data=data.frame(shade_cent=s, water_cent=-1:1))
    for(i in 1:20) lines(-1:1,  mu[i, ], col=col.alpha("black", .3), lwd=ifelse(i==1, 3, 1))
  }
}
