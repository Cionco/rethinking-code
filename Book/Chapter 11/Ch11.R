library(rethinking)
data(chimpanzees)
d <- chimpanzees


d$treatment <- 1 + d$prosoc_left + 2 * d$condition

xtabs(~ treatment + prosoc_left + condition, d)

m11.1 <- quap(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a,
    a ~ dnorm(0, 1.5)
  ), data=d
)


set.seed(1999)
prior <- extract.prior(m11.1, n=1e4)

p <- inv_logit(prior$a)
dens(p, adj=.1)

m11.2 <- quap(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a + b[treatment],
    a ~ dnorm(0, 1.5),
    b[treatment] ~ dnorm(0, 0.5)
  ), data=d
)

prior <- extract.prior(m11.2, n=1e4)
p <- sapply(1:4, function(k) inv_logit(prior$a + prior$b[, k]))

dens(abs(p[,1]-p[,2]), adj=.1)




m11.3 <- quap(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a + b[treatment],
    a ~ dnorm(0, 1.5),
    b[treatment] ~ dnorm(0, 0.5)
  ), data=d
)

prior <- extract.prior(m11.3, n=1e4)
p <- sapply(1:4, function(k) inv_logit(prior$a + prior$b[, k]))

mean(abs(p[,1]-p[,2]))



#trimmed data list
dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = as.integer(d$treatment)
)

m11.4 <- ulam(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a[actor] + b[treatment],
    a[actor] ~ dnorm(0, 1.5),
    b[treatment] ~ dnorm(0, 0.5)
  ), data=dat_list, chains=4, log_lik=TRUE
)
precis(m11.4, depth=2)

post <- extract.samples(m11.4)
p_left <- inv_logit(post$a)
plot(precis(as.data.frame(p_left)), xlim=c(0, 1))

labs <- c("R/N", "L/N", "R/P", "L/P")
plot(precis(m11.4, depth=2, pars="b"), labels=labs)

diffs <- list(
  db13 = post$b[,1] - post$b[,3],
  db24 = post$b[,2] - post$b[,4]
)
plot(precis(diffs))











pl <- by(d$pulled_left, list(d$actor, d$treatment), mean)
pl[1,]

plot( NULL,xlim=c(1,28),ylim=c(0,1),xlab="",
      ylab="proportion leftlever",xaxt="n",yaxt="n")
axis( 2,at=c(0,0.5,1),labels=c(0,0.5,1))
abline( h=0.5,lty=2)
for (jin1:7)abline(v=(j-1)*4+4.5,lwd=0.5)
for (jin1:7)text((j-1)*4+2.5,1.1,concat("actor",j),xpd=TRUE)
for (jin(1:7)[-2]){
  lines( (j-1)*4+c(1,3),pl[j,c(1,3)],lwd=2,col=rangi2)
  lines( (j-1)*4+c(2,4),pl[j,c(2,4)],lwd=2,col=rangi2)
}
points( 1:28,t(pl),pch=16,col="white",cex=1.7)
points( 1:28,t(pl),pch=c(1,1,16,16),col=rangi2,lwd=2)
yoff <-0.01
text( 1,pl[1,1]-yoff,"R/N",pos=1,cex=0.8)
text( 2,pl[1,2]+yoff,"L/N",pos=3,cex=0.8)
text( 3,pl[1,3]-yoff,"R/P",pos=1,cex=0.8)
text( 4,pl[1,4]+yoff,"L/P",pos=3,cex=0.8)
mtext( "observedproportions\n")









dat <- list(actor=rep(1:7, each=4), treatment=rep(1:4, times=7))
p_post <- link(m11.4, data=dat)
p_mu <- apply(p_post, 2, mean)
p_ci <- apply(p_post, 2, PI)




plot( NULL,xlim=c(1,28),ylim=c(0,1),xlab="",
      ylab="proportion leftlever",xaxt="n",yaxt="n")

axis(2, at=c(0, 0.5, 1), labels=c(0, 0.5, 1))

abline(h=0.5, lty=2)
for (j in 1:7) abline(v = (j - 1) * 4 + 4.5, lwd=0.5)
for (j in 1:7) text((j - 1) * 4 + 2.5, 1, concat("actor",j), xpd=TRUE)
for (j in (1:7)[-2]) {
  index_vectors <- cbind((j - 1) * 4 + c(1, 3), (j - 1) * 4 + c(2, 4))

  lines(index_vectors[, 1], p_mu[index_vectors[, 1]], lwd=2, col="black")
  lines(index_vectors[, 2], p_mu[index_vectors[, 2]], lwd=2, col="black")
}
for(i in 1:28) {
  lines(rep(i, 2), p_ci[,i], lwd=2)
}

points( 1:28, p_mu, pch=16, col="white", cex=c(1, 1, 1.1, 1.1))
points( 1:28, p_mu, pch=c(1,1,16,16), col="black", lwd=1, cex=c(1, 1, 0.8, 0.8))

mtext("posterior predictions")








