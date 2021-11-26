library(rethinking)

data(WaffleDivorce)
d <- WaffleDivorce

d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)


# Relationship Divorce - Median Age
m5.1 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bA * A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)

set.seed(10)
prior <- extract.prior(m5.1)
mu <- link(m5.1, post=prior, data=list(A=c(-2, 2)))
plot(NULL, xlim=c(-2, 2), ylim=c(-2, 2))
for(i in 1:50) lines(c(-2, 2), mu[i,], col=col.alpha("black", 0.4))

seq <- seq(from=-3, to=3.2, length.out=30)
mu <- link(m5.1, data=list(A=seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(D ~ A, data=d, col=rangi2)
lines(seq, mu.mean, lwd=2)
shade(mu.PI, seq)


# Relationship Divorce - Marriage Rate
m5.2 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)

prior <- extract.prior(m5.2)
mu <- link(m5.2, post=prior, data=list(M=c(-2, 2)))
plot(NULL, xlim=c(-2, 2), ylim=c(-2, 2))
for(i in 1:50) lines(c(-2, 2), mu[i,], col=col.alpha("black", 0.4))

mu <- link(m5.2, data=list(M=seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(D ~ M, data=d, col=rangi2)
lines(seq, mu.mean, lwd=2)
shade(mu.PI, seq)


# Multiple Regression
m5.3 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M + bA * A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)

plot(coeftab(m5.1, m5.2, m5.3), par=c("bA", "bM"))


# Marriage Rate ~ Marriage Age

m5.4 <- quap(
  alist(
    M ~ dnorm(mu, sigma),
    mu <- a + bAM * A,
    a ~ dnorm(0, 0.2),
    bAM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data=d
)

mu <- link(m5.4)
mu.mean <- apply(mu, 2, mean)
mu.resid <- d$M - mu.mean


plot(mu.resid, d$D)



mu <- link(m5.3)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

D.sim <- sim(m5.3, n=1e4)
D.PI <- apply(D.sim, 2, PI)

plot(mu.mean ~ d$D, col=rangi2, ylim=range(mu.PI), xlab="Observed Divorce", ylab="Predicted Divorce")
abline(a=0, b=1, lty=2)
for(i in 1:nrow(d)) lines(rep(d$D[i], 2), mu.PI[,i], col=rangi2)
identify(x=d$D, y=mu.mean, labels=d$Loc)



m5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M + bA * A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1),
    
    ##A -> M
    M ~ dnorm(mu_M, sigma_M),
    mu_M <- aM + bAM * A,
    aM ~ dnorm(0, 0.2),
    bAM ~ dnorm(0, 0.5),
    sigma_M ~ dexp(1)
  ), data=d
)

A_seq <- seq(-2, 2, length.out=30)

sim_dat <- data.frame(A=A_seq)
s <- sim(m5.3_A, data=sim_dat, vars=c("M", "D"))

plot(sim_dat$A, colMeans(s$D), ylim=c(-2, 2), type="l", xlab="manipulated A", ylab="counterfactual D")
shade(apply(s$D, 2, PI), sim_dat$A)
mtext("Total counterfactual effect of A on D")


sim_dat <- data.frame(M=seq(-2, 2, length.out=30), A=0)
s <- sim(m5.3_A, data=sim_dat, vars="D")

plot(sim_dat$M, colMeans(s), ylim=c(-2, 2), type="l", xlab="manipulated M", ylab="counterfactual D")
shade(apply(s, 2, PI), sim_dat$M)
mtext("Total counterfactual effect of M on D")
