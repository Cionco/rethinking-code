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

