library(stargazer)
library(lmtest)
library(sandwich)
library(car)
library(margins)


load("Data20745516.rda")

PartA<-na.omit(Hours)

# Creating Duumy Variable "haveKids"
PartA$havekids <- ifelse(Hours$youngkids > 0, 1, 0)
PartA$havekids <- ifelse((Hours$oldkids > 0), 1, PartA$havekids)

# The 3 models are:

model_1 <- lm(hours~youngkids+I(youngkids^2)+I(havekids*hhours^2)+havekids+log(age)+experience+hhours+log(hwage),data=PartA)
model_2 <- lm(log(hours)~youngkids+I(havekids*hhours)+havekids+experience+age+hhours+I(hhours*hwage^2)+log(hwage),data=PartA,subset=(hours>0))
model_3 <- lm(hours~youngkids+I(havekids*hhours^2)+havekids+log(age)+experience+hhours+log(hwage),data=PartA)

# from Proj2 file
printEqu <- function(obj,digits=4, form=NULL, maxpl=5, adjrsq=FALSE, 
                     robsd=NULL, stars=FALSE, dist=c("t","n"), label=NULL,
                     omit=NULL, ...)
{
  dist <- match.arg(dist)
  if (is.null(label))
  {
    typeEQ <- "equation*"
    label=""
  } else {
    typeEQ <- "equation"
    label <- paste("\\label{",label,"}", sep="")
  }
  if (!is.null(omit))
  {
    omit <- omit[omit != "(Intercept)"]
    chk <- lapply(omit, function(o) grep(o, names(obj$coef)))
    chk <- sapply(1:length(chk), function(i) length(chk[[i]])==0)
    omit <- omit[!chk]
    if (length(omit)==0)
      omit <- NULL
  }
  if (is.null(form))
  {                
    cat("\\begin{",typeEQ, "}", label,"\n", sep="")
    cat("\\begin{split}\n")
    ncoef <- names(coef(obj))
    if (!is.null(omit))
    {
      omit <- do.call("c", lapply(omit, function(o) grep(o, ncoef)))
      omit <- -unique(omit)
      nomit <- length(omit)
      add <- paste("\\mbox{ (+ ", nomit, " omitted terms)}", sep="")
    } else {
      add <- ""
      omit <- 1:length(ncoef)
    }                
    Intercept <- attr(obj$terms, "intercept")
    b <- formatC(abs(obj$coef), digits=digits, ...)[omit]
    if (!is.null(robsd))
      snum <- robsd[omit]
    else
      snum <- summary(obj)$coef[omit,2]
    s <- formatC(snum, digits=digits, ...)
    ncoef <- ncoef[omit]
    if (stars)
    {
      ttest <- coef(obj)[omit]/snum
      if (dist=="t")
        pv <- 2*pt(-abs(ttest), obj$df)
      else
        pv <- 2*pnorm(-abs(ttest))
      sym <- symnum(pv, cutpoints=c(0,.01,.05,.1,1),
                    symbols=c("^{***}","^{**}","^*"," "))
      symmess <- "\\\\& ^*\\text{pv}<0.1\\mbox{; }^{**}\\text{pv}<0.05\\mbox{; }^{***}\\text{pv}<0.01"
    } else {
      sym <- rep("", length(coef(obj)[omit]))
      symmess <- ""
    }                                   
    ny <- rownames(attr(obj$terms, "factors"))[1]                     
    ny <- paste("\\widehat{",ny,"}",sep="")
    cat(ny,"&=")
    if (obj$coef[1] < 0)
      cat("\\underset{(",s[1],")",sym[[1]],"}{-",b[1],"}", sep="")
    else
      cat("\\underset{(",s[1],")",sym[[1]],"}{",b[1],"}")                    
    if (Intercept==0)
      cat("~",ncoef[1], sep="")
    j <- 1
    for (i in 2:length(b))
    {
      if (j>maxpl)
      {
        j <- 1
        cat("\\\\&\\quad\n")
      }
      if ((obj$coef[omit])[i] < 0)
        cat("~-~")
      else
        cat("~+~")
      cat("\\underset{(",s[i],")",sym[[i]],"}{",b[i],"}~")
      cat(ncoef[i])
      j <- j+1    
    }
    cat(add)
    n <- length(obj$residuals)
    cat("\\\\ &\\quad n=", n, ",~~R^2=", round(summary(obj)$r.squared,5))
    cat(", SSR=", round(sum(obj$resid^2),5))
    if (adjrsq)
      cat(", \\bar{R}^2=", round(summary(obj)$adj,5))
    if (!is.null(robsd))
      cat("\\mbox{ (Robust S-E)}")
    cat(symmess)
    cat("\n\\end{split}\n")
    cat("\\end{", typeEQ, "}\n", sep="")
  } else {
    t <- terms(form)
    y <- rownames(attr(t, "factors"))[1]                     
    x <- colnames(attr(t, "factors"))
    cat("\\begin{",typeEQ, "}", label,"\n\\begin{split}\n",
        y,"&=\\beta_0", sep="")
    j  <-  1
    for (i in 1:length(x))
    {
      if (j>maxpl)
      {
        j <- 1
        cat("\\\\&\n")
      }
      cat("+\\beta_",i,x[i],sep="")
      j <- j+1
    }
    cat("+u\n\\end{split}\n\\end{",typeEQ,"}", sep="")
  }
}

printEqu(model_1, stars=TRUE, maxpl=3,adjrsq=TRUE, digits=5)
printEqu(model_2, stars=TRUE, maxpl=3,adjrsq=TRUE, digits=5)
printEqu(model_3, stars=TRUE, maxpl=3,adjrsq=TRUE, digits=5)

#We first estimate the models and test the homoscedasticity using the short White test:

yhat <- fitted(model_1)
yhat2 <- yhat^2
bptest(model_1, ~yhat+yhat2)

yhat <- fitted(model_2)
yhat2 <- yhat^2
bptest(model_2, ~yhat+yhat2)

yhat <- fitted(model_3)
yhat2 <- yhat^2
bptest(model_3, ~yhat+yhat2)

# Indirect t-test
coeftest(model_1, vcov.=vcovHC, df=Inf)

# J-test between Model 1 and 2
jtest(model_1, model_2,vcov=vcovHC, df=Inf, data=PartA)

# Checking if model properly specified
resettest(model_1, power=2:3)

bp <- bptest(model_1)
adj <- if(bp$p.value<.05) vcovHC else NULL
coeftest(model_1, vcov.=adj, df=Inf)

#Cook's Distance
plot(model_1,4)

# Removing outliers
fit <- lm(hours~youngkids+I(youngkids^2)+I(havekids*hhours)+havekids+
            age+experience+hhours+hwage,data=PartA,
          subset = !(rownames(PartA) %in% c(37,126,403,598,734,292))) 
plot(fit,4)

# impact of dropping these three observations:

stargazer(model_1, fit, type='text')

confint(fit)

# APE of youngkids on hours
b <- coef(fit)
m <- mean(PartA$youngkids)
APE <- b[2]+2*b[3]*m

APE
# 95% Confidence of the APE
v <- vcov(fit)
s <- sqrt(v[2,2]+4*m^2*v[3,3]+4*m*v[2,3])

crit <- qt(.975, fit$df)
c(APE-s*crit, APE+s*crit)

#  Does having children less than 6 have the same effect as having
#  children between ages 6 and 18?
linearHypothesis(fit,c("youngkids","I(youngkids^2)"))


## Women with kids VS Women without kids

fit <- lm(hours~youngkids+I(youngkids^2)+I(havekids*hhours)+havekids+
            age+experience+hhours+hwage,data=PartA,
          subset = !(rownames(PartA) %in% c(37,126,403,598,734,292)))

new<-PartA[!(rownames(PartA) %in% c(37,126,403,598,734,292)),]

r <- range(new$hhours)

hhours <- r[1]:r[2]
newd1 <- data.frame(hhours=hhours, youngkids= mean(new$youngkids),
                    havekids=1,  age= mean(new$age),
                    experience= mean(new$experience),hwage=mean(new$hwage))
newd2 <- data.frame(hhours=hhours, youngkids= mean(new$youngkids),
                    havekids=0,  age= mean(new$age),
                    experience= mean(new$experience),hwage=mean(new$hwage))

pr1 <- predict(fit, newdata=newd1)
pr2 <- predict(fit, newdata=newd2)
pred <- cbind(pr1,pr2)

col <- (new$havekids==1) +
  2*(new$havekids==0)
pch <- 21*(new$havekids==1) +
  22*(new$havekids==0) 
plot(hours~hhours, col=col, bg=col, pch=pch, data=new)
matplot(hhours, pred, col=1:2, lty=1:2, lwd=2, type='l',
        main="Predicted working hours of women with and without kids in 1975",
        ylab="hours",xlab="husband's working hours", add=TRUE)
legend(x=3500, y=2900,legend=c("women with kids", "women without kids"),col=1:2,
       lty=1:2, lwd=2,pt.bg=1:4, pch=21:24, cex=.8)
grid()
