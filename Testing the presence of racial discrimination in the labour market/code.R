library(stargazer)
library(lmtest)
library(sandwich)
library(car)
library(margins)

load("reg_data.rda")

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

PartB<-Names

# Changinging categorical data to binary 1 and 0.

PartB$gender <- as.numeric(PartB$gender == "male")
PartB$call <- as.numeric(PartB$call == "yes")
PartB$ethnicity <- as.numeric(PartB$ethnicity == "afam")
PartB$email <- as.numeric(PartB$email == "yes")
PartB$military <- as.numeric(PartB$military == "yes")
PartB$volunteer <- as.numeric(PartB$volunteer == "yes")
PartB$equal <- as.numeric(PartB$equal == "yes")
PartB$city <- as.numeric(PartB$city == "chicago") # Chicago is 1 and boston is 0


# Estimating the models
fit_1 <- lm (call~ ethnicity + city + I(city * ethnicity) + experience +
               I(experience^2),data = PartB)

fit_2 <- lm (call ~ ethnicity+ gender + military + ethnicity * military
             + (gender * military) + (ethnicity * gender)+
               (ethnicity* gender *military) + experience +
               I(experience^2),data=PartB )

fit_3 <- lm(call ~(email*ethnicity)*(experience+I(experience^2)),data=PartB)
fit_4 <- lm(call ~ volunteer + ethnicity + equal+ I(ethnicity * equal) +
              I (volunteer * equal) + I (ethnicity* volunteer)+ experience +
              I(experience^2),data=PartB)

# QUESTION 1

# Printing coefficient table with robust s.e. and test 

coeftest(fit_1, vcov.=vcovHC)

# Testing for African-American sounding names in Chicago VS Caucasian sounding names in Chicago: 
# Using the robust F test:

linearHypothesis(fit_1, "ethnicity + I(city * ethnicity)=0", white.adjust=TRUE, test="Chisq")

# Testing for African-American sounding names in Boston VS Caucasian sounding names in Boston: 

# Using the robust F test:

# ethnicity =0
linearHypothesis(fit_1, "ethnicity = 0", white.adjust=TRUE, test="Chisq")

# Testing for African-American sounding names in Chicago VS African-American sounding names in Boston: 

# Using the robust F test:

# "city + I(city * ethnicity) = 0"
linearHypothesis(fit_1, "city + I(city * ethnicity) = 0",
                 white.adjust=TRUE, test="Chisq")

# Computing robust confidence intervals:

coefci(fit_1, level=.95, vcov.=vcovHC, df=Inf)

## QUESTION 2

# Testing Condition 1:

# Using the robust F test:

# "ethnicity + ethnicity:military+ethnicity:gender+ethnicity:gender:military= 0"
linearHypothesis(fit_2, "ethnicity + ethnicity:military+ethnicity:gender+
                 ethnicity:gender:military= 0",
                 white.adjust=TRUE,vcov.=vcovHC, test="Chisq")

# Testing Condition 2:

# Using the robust F test:

# "ethnicity+ethnicity:military =0"
linearHypothesis(fit_2, "ethnicity+ethnicity:military =0",
                 white.adjust=TRUE,vcov.=vcovHC, test="Chisq")

# Testing Condition 3:

# Using the robust F test:

# "ethnicity+ ethnicity:gender=0"
linearHypothesis(fit_2, "ethnicity+ ethnicity:gender=0",
                 white.adjust=TRUE,vcov.=vcovHC, test="Chisq")

# Testing Condition 4:

# Using the robust F test:

# "ethnicity=0"
linearHypothesis(fit_2, "ethnicity=0",
                 white.adjust=TRUE,vcov.=vcovHC, test="Chisq")

# Testing Condition 5:

# Using the robust F test:
# "gender +gender:military+ethnicity:gender+ethnicity:gender:military=0"
linearHypothesis(fit_2, "gender+gender:military+ethnicity:gender+
                 ethnicity:gender:military=0",
                 white.adjust=TRUE,vcov.=vcovHC, test="Chisq")

# Testing Condition 6:

# Using the robust F test:

# "gender + ethnicity:gender= 0"
linearHypothesis(fit_2, "gender + ethnicity:gender= 0",
                 white.adjust=TRUE,vcov.=vcovHC, test="Chisq")

# Computing robust confidence intervals:

coefci(fit_2, level=.95, vcov.=vcovHC, df=Inf)

# QUESTION 3

# Printing coefficient table with robust s.e. and test 

coeftest(fit_3, vcov.=vcovHC)

# Testing Condition 1:

# Using the robust F test:

# "ethnicity+email:ethnicity+ethnicity:experience+ethnicity:experience^2+
#  email:ethnicity:experience+email:ethnicity:I(experience^2)= 0"

linearHypothesis(fit_3, "ethnicity+email:ethnicity+
                 ethnicity:experience+ethnicity:I(experience^2)+
                 email:ethnicity:experience+
                 email:ethnicity:I(experience^2)= 0",
                 white.adjust=TRUE, test="Chisq")

# Testing Condition 2:

# Using the robust F test:

# "ethnicity+ethnicity:experience+ethnicity:I(experience^2)=0"

linearHypothesis(fit_3, "ethnicity+ethnicity:experience+
                 ethnicity:I(experience^2)=0",white.adjust=TRUE, test="Chisq")

# Computing robust confidence intervals:

coefci(fit_3, level=.95, vcov.=vcovHC, df=Inf)


## QUESTION 4

# Printing coefficient table with robust s.e. and test 

coeftest(fit_4, vcov.=vcovHC)

# Testing Condition 1:

linearHypothesis(fit_4, "ethnicity+I(ethnicity * equal)+I(ethnicity * volunteer)= 0", white.adjust=TRUE, test="Chisq")

# Testing Condition 2:

linearHypothesis(fit_4, "ethnicity + I(ethnicity * equal) = 0", white.adjust=TRUE, test="Chisq")

# Testing Condition 3:

linearHypothesis(fit_4, "ethnicity+I(ethnicity * volunteer)=0", white.adjust=TRUE, test="Chisq")

# Testing Condition 4:

linearHypothesis(fit_4, "ethnicity = 0", white.adjust=TRUE, test="Chisq")

# Computing robust confidence intervals:

coefci(fit_4, level=.95, vcov.=vcovHC, df=Inf)


library(margins)
a<-margins_summary(fit_3, vcov=vcovHC(fit_3))

# The Average Partial effect of experience on call is:

a$AME[3]

# The confidence interval for this is:

c(a$lower[3],a$upper[3])

# The Average Partial effect of email on call is:

a$AME[1]

# The confidence interval for this is:

c(a$lower[1],a$upper[1])

# The Average Partial effect of ethnicity on call is:

a$AME[2]

# The confidence interval for this is:

c(a$lower[2],a$upper[2])

# Graph of different groups
newd <- data.frame(email=1, experience=seq(0,30),ethnicity=1)
pr <- predict(fit_3, newdata=newd)
newd <- data.frame(email=1, experience=seq(0,30),ethnicity=0)
pr2 <- predict(fit_3, newdata=newd)
newd <- data.frame(email=0, experience=seq(0,30),ethnicity=1)
pr3 <- predict(fit_3, newdata=newd)
newd <- data.frame(email=0, experience=seq(0,30),ethnicity=0)
pr4 <- predict(fit_3, newdata=newd)
plot(call~experience, data=PartB,
     main="Probability of calling back versus experience",ylab="call",
     xlab="experience",ylim=c(0,1.7))
lines(newd[,2], pr, col=2, lwd=2)
lines(newd[,2], pr2, col=4, lwd=2)
lines(newd[,2], pr3, col="green", lwd=2)
lines(newd[,2], pr4, col="yellow", lwd=2)
legend("topright",c("African-American with email","Caucasian with email",
                    "African-American with no email","Caucasian with no email"),
       col=c(2,4,"green","yellow"), lty=2,lwd=2)

# Confidence Interval
newd <- data.frame(email=1, experience=seq(0,30),ethnicity=1)
predict(fit_3, newdata=newd,interval = "confidence")
newd <- data.frame(email=1, experience=seq(0,30),ethnicity=0)
predict(fit_3, newdata=newd,interval = "confidence")
newd <- data.frame(email=0, experience=seq(0,30),ethnicity=1)
predict(fit_3, newdata=newd,interval = "confidence")
newd <- data.frame(email=0, experience=seq(0,30),ethnicity=0)
predict(fit_3, newdata=newd,interval = "confidence")