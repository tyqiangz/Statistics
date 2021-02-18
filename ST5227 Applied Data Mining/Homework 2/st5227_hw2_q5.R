xy = read.table('DMleukemia1.dat')
xy = data.matrix(xy)

y = xy[,251]
y = (y>0) - (y<=0)   # redefine the value of the classes 

n = length(y)
X = xy[,1:250]

Xmean = apply(X, 2, mean)
Xstd = apply(X, 2, sd)
ymean = mean(y);

X = scale(X)
y = y-ymean

I = diag(250);

cv = matrix(0, 100, 1) 
lambda = ((1:100)/100)^5/10
for (i in 1:100)
{     
  lambdai = lambda[i] # value of lambda
  cvi = 0
  for (j in 1:n)
  {
    bR = solve( t(X[-j,]) %*% X[-j,] + lambdai*I) %*% t(X[-j,]) %*% y[-j]; # ridge estimator of beta
    yje = X[j,] %*% bR  # fitted values
    cvi = cvi + ((yje-y[j])^2)
  }
  cv[i] = cvi/n
}
plot(log(lambda), log(cv)) 

BestLambda = lambda[which(cv==min(cv))]

BestLambda

Ridge = solve( t(X) %*% X + BestLambda*I) %*% t(X) %*% y; # ridge estimator of beta


######## model validation for a new data set

xyNEW = read.table('DMleukemia2.dat')
xyNEW = data.matrix(xyNEW)
yNEW = xyNEW[,251]
yNEW = (yNEW>0) - (yNEW<=0)

n = length(yNEW)
XNEW = xyNEW[,1:250]
XNEW  = XNEW - matrix(rep(Xmean, each=n), n, 250)
XNEW = XNEW/matrix(rep(Xstd, each=n), n, 250)


ypredicted = XNEW%*%Ridge + ymean

errorRidge = mean((yNEW-ypredicted)^2)
errorRidge

# classification error

yclassified = (ypredicted>0)- (ypredicted<=0)
ClassficationError = mean(yNEW != yclassified)

ClassficationError

plot(ypredicted, yNEW, pch = 20, ylab='classes', xlab='predicted values')
points(ypredicted, yclassified, col='red', pch=2) 
legend(-1, 0.9, pch=20, "true")
legend(-1, 0.4, pch=2, "classified", col='red')

############################ path

lambda = (1:100)^4
nL = length(lambda)
BETA = matrix(0, nL, 250)
for (i in 1:nL)
{     
  lambdai = lambda[i] 
  BETA[i,] = solve( t(X) %*% X + lambdai*I) %*% t(X) %*% y; 
}

t = apply(BETA^2, 1, sum)
matplot(t, BETA, type='l')
