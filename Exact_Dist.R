#See http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2476536
#Andy Wheeler, apwheele@gmail.com

#functions used in ExactProb - note these do no error checking 
#for observations in zero probability bins - so can return non-sensical results
#Minimalist chi square, default equal probability for each bin
chiStat <- function(v,p=rep(1/length(v),length(v))){
  n <- sum(v)
  e <- n*p
  r <- (v-e)^2
  c <- ifelse(e>0,r/e,0) #this is for zero prob bins for alternate process
  return(sum(c))
}
#Kuipers V - probably stole this from the circular package at some point
KuiperTest <- function(v,p=rep(1/length(v),length(v))){
  n <- sum(v)
  u <- cumsum(p) 
  s <- v/n
  e <- cumsum(s)        #ecdf
  Dp <- max(e - u)
  Dm <- max(s)
  sq_n <- sqrt(n)
  V <- (Dp + Dm)* (sq_n + 0.155 + 0.24/sq_n)
  return(V) 
}
#Kolmogrov smirnov test
KSTest <- function(v,p=rep(1/length(v),length(v))){
  n <- sum(v)
  u <- cumsum(p)
  s <- v/n
  e <- cumsum(s)        #ecdf
  Dp <- max(e - u)
  Dm <- max(s)
  D <- max(c(Dp,Dm))
  return(D) 
}
#likelihood ratio G test
GTest <- function(v,p=rep(1/length(v),length(v))){
  e <- sum(v)*p
  r <- ifelse(v>0&e>0,log(v/e),0)
  g <- 2*sum(v*r)
  return(g)
}
#multinomial prob based on set of probabilities, defaults to equal probabilities
exactMult <- function(v,p=rep(1/length(v),length(v))){
    n <- factorial(sum(v))
    d <- prod(factorial(v))
    p <- prod(p^v)
    return( (n/d)*p )
}
#This generates all the permutations given n number of balls in m bins and 
#then calculates the exact probability according to the multinomial 
#distribution and the CDF of the specified stat
exactProb <- function(n,m,p=rep(1/m,m),type="Chi"){
  library(partitions)
  AllDat <- t(compositions(n,m))
  ExactProb <- apply(AllDat,1,exactMult,p=p)
  if (type == "Chi"){
    Stat <- apply(AllDat,1,chiStat,p=p) }
  else if (type == "V") {
    Stat <- apply(AllDat,1,KuiperTest,p=p) }
  else if (type == "G") {
    Stat <- apply(AllDat,1,GTest,p=p) }
  else if (type == "KS") {
    Stat <- apply(AllDat,1,KSTest,p=p) }
  #order according to stat
  MyData <- data.frame(as.matrix(AllDat),ExactProb, Stat)[order(Stat),]
  #MyData$cumprob <- cumsum(MyData$ExactProb)
  return(MyData)
}
#My wrapping all up in a global function to return items in list
#given the initial data
#Options are likelihood ratio "G" test, Kuiper's "V" for circular data
#"KS" for Kolmogrov smirnov test
#and "Chi" for chi-square, I found G for my paper was always more powerful when 
#total n>8 for 7 bins, and was equal in power to Chi-square below
#Kuipers V was more powerful when using circular data in some circumstances

SmallSampTest <- function(d,p=rep(1/length(d),length(d)),type="G"){
  n <- sum(d)
  m <- length(d)
  cdf <- exactProb(n=n,m=m,p=p,type=type)  #generate exact probability
  if (type == "Chi"){
    Samp <- chiStat(d,p) }
  else if (type == "V") {
    Samp <- KuiperTest(d,p)  }
  else if (type == "G") {
    Samp <- GTest(d,p)  }
  else if (type == "KS") {
    Samp <- KSTest(d,p)  }  
  #p-value to the right of the test statistic, aggregating first
  Agg <- aggregate(x=cdf[,'ExactProb'],by=list(cdf$Stat),FUN=sum)
  names(Agg) <- c("Stat","ExactProb")
  Agg$cumprob <- cumsum(Agg$ExactProb)
  pvalue <- sum(Agg[Agg[,'Stat'] >= Samp,'ExactProb'])
  #return object
  t <- list(cdf,p,d,Samp,pvalue,Agg)
  names(t) <- c("CDF","probabilities","data",type,"p-value","Aggregate Statistics")
  class(t) <- "SmallSampleTest"
  return(t)
}
#calculating the power of an alternative test
PowAlt <- function(SST,p_alt,a=.05){
  x <- merge(SST$CDF,SST$`Aggregate Statistics`,by="Stat")
  x$AltProb <- apply(as.matrix(x[,2:(1+length(p_alt))]),1,exactMult,p=p_alt)
  power <- sum(x[x[,'cumprob'] > (1-a),'AltProb']) #power of alt
  r <- list(x,power,p_alt,SST$probabilities,a)
  names(r) <- c("permutations","power","Alternative","null","alpha")
  names(r$permutations)[1] <- names(SST[4])
  class(r) <- "PowerSmallSamp"
  return(r)
}
#also see http://stats.stackexchange.com/a/125150/1036
#for an example of calculating the power under a particular alternative

#print function for class 
print.SmallSampleTest <- function(x,...){
  cat("Small Sample Test Object \n")
  cat(paste0("Test Type is ",names(x[4])," \n"))
  cat(paste0("Statistic is ",x[4]," \n"))
  cat("p-value is: ",x$'p-value'," \n")
  cat("Data are: ",paste(x$data),"\n")
  cat("Null probabilities are: ",formatC(x$probabilities,digits=2),"\n")
  cat("Total permutations are: ",length(x$CDF[,1])," \n")
}

print.PowerSmallSamp <- function(x){
  cat("Power for Small Sample Test \n")
  cat("Test statistic is:",names(x$permutations)[1]," \n")
  cat("Power is:",x$power," \n")
  cat("Null is:",formatC(x$null,digits=2)," \n")
  cat("Alt is:",formatC(x$Alternative,digits=2)," \n")
  cat("Alpha is:",x$alpha," \n")
  b <- length(names(x$permutations))-5
  cat("Number of Bins:",b," \n")
  s <- sum(x$permutations[1,2:(1+b)])
  cat("Number of Observations:",s," \n")
}
 
#TODO
#error checking for power function
#pass your own function
#example power analysis
#R package

##now with an example dataset, http://stats.stackexchange.com/q/133377/1036
#d <- c(3,2,1,2,1,2,6) #format N observations in M bins
#p <- rep(1/7,7) #equiprobable across days
#t <- SmallSampTest(d=d,p=p,type="Chi")
#t
##using the asymptotic distribution you come to the same conclusion
#chisq.test(d) #p-value slightly lower
##using simulation p-value is closer to exact function
#chisq.test(d,simulate.p.value = TRUE, B = 9999)

##likelihood ratio G statistic
#gT <- SmallSampTest(d=d,p=p,type="G")
#gT
#
##using the asymptotic distribution you come to the same conclusion
#chisq.test(d)
#
##power of the test under alternative hypothesis
#x <- c(1,0,4)
#gV <- SmallSampTest(d=x,type="G")
#p_alt <- c(0.1,0.1,0.8)
#PowAlt(SST=gV,p_alt=p_alt)
#
#SmallSampTest(d=x,type="Chi")
#chisq.test(x)
#
##Example with KS test
#mC <- c(1,2,3)
#SmallSampTest(d=mC,type="KS")
##same D value as ks test with expanded data
##because CDF of uniform is the same as equal prob bins
#x <- c(1,2,2,3,3,3)
#ks.test(x,"punif",min=1,max=3)


