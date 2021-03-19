#scatter plot
attach(n)
plot(pcc, nd, main="scatplot", xlab="pcc", ylab="network density", pch=19)
abline(lm(nd~pcc), col="red") # regression line (y~x)