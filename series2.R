# This function computes the odds P(nochange)/P(change) in two binomial sequences
# n1,n2 are the sample sizes of the sequences
# x1,x2 is the counting of 1s in the sequences
# The problem is solved with a Bayesian BDEu score with equivalent sample size s

computeoddsc <- function(n1,n2,x1,x2,s) {
  
  lgs4 <- lgamma(s/4)
  lgs2 <- lgamma(s/2)
  
change <-(lgamma(s/2) - lgamma(n1+s/2)     + 
            lgamma(x1+s/4 ) -  lgs4  +
            lgamma((n1-x1)+s/4 ) -  lgs4 +
          lgamma(s/2) - lgamma(n2+s/2)     + 
            lgamma(x2+s/4 ) -  lgs4  +
            lgamma((n2-x2)+s/4 ) -  lgs4
)

nochange <-  (lgamma(s) - lgamma(n1+n2+s) +
                lgamma(x1+x2 +s/2)  - lgs2 +
                lgamma(n1+n2-x1-x2+s/2) - lgs2
)


relative <- nochange-change


odd <- exp(relative) 
return(odd)

}

#simulate a series of binomial values with with t a vector
# of the different probability values and n the size of the sample for 
# each probability value (the same for all)

simulate <- function(n,t) {
x <- c()
  
for (r in t) {

x1<- rbinom(n,1,r)
x<-append(x,x1)

}
return(x)
}


# Function that estimates probabilities from a string x
# It returns a list with the estimations, the sample sizes, and the forgotten samples
# It contains a window of active values. It has a parameter n 
# for making the test: it compares the n last values with the n first
# values of the windos (it must have a size greater or equal to 2*n)
# if there are significant differences it removes the first n values
# of the window and repeats the test.

estimate1 <- function(x,n) {
  
y <- vector(,length(x))  
s   <- vector(,length(x))  
ro <-  vector(,length(x))  
l<-0

k<-1

for(i in 1:length(x)){
  
  if (i-k>=2*n-1) {
    test<- 1
    l<- 0
    
    while(test==1) {
      j1 = k+n-1
      j2 = i-n+1
      
      x1 = sum(x[k:j1])
      x2 = sum(x[j2:i])
      
      odd <- computeoddsc(n,n,x1,x2,4)
      
      if (odd <0.5) {k<-k+n
      l<- l+n
      if (i-k<2*n-1) {test<-0} 
      }
      else {
        test<-0
      }
    }
    
  }
  
  y[i]<-(sum(x[k:i])+1)/(i-k+1+2)
  s[i] <-  i-k+1
  ro[i] <- l
  
}

return(list(y,s,ro))

}

   

# Function that estimates probabilities from a string x
# It returns a list with the estimations, the sample sizes,  the forgotten samples, and the forgetting coefficients
# It makes statistical tests (goodness of fit) of the n1 last samples with respect to the
# estimations done on step i-n1. If test is meaningfull, then it multiplies the past samples by 
# a coefficient rho


estimate2 <- function(x,n1){

y <- vector(,length(x))  
s   <- vector(,length(x))  
ro <-  vector(,length(x))  
fg <- vector(,length(x)) 
r <- vector(,length(x))  

r[1]=x[1]
ro[1] = 1
s[1] = 1
y[1] = (r[1]+1)/(s[1]+2.0)
fg[1] = 0

for(i in 2:length(x)){

  if (i>n1) {
   
    pvalue <- binom.test(sum(x[(i-n1+1):i]),n1,y[i-n1])$p.value
  if(pvalue>0.05) {ro[i]<-1
  }
  else {ro[i] <- pvalue^{1/100}}
  }
  else {ro[i]<-1}
  r[i] <- r[i-1]*ro[i]+x[i]
  fg[i] <- s[i-1] * (1-ro[i])
  s[i]<-  s[i-1]*ro[i]+1
  y[i]<- (r[i]+1)/(s[i]+2.0)
  
}

return(list(y,s,fg,ro))

}

# Function that estimates the probabilities with a constant forgeting factor ro

estimate3 <- function(x,ro) {
  
  y <- vector(,length(x))  
  s   <- vector(,length(x))  
  rov <- rep(ro,length(x))
  fg <- vector(,length(x)) 
  r <- vector(,length(x))  
  
r[1] = x[1]
s[1] = 1
fg[1]=0
for(i in 2:length(x)){
  r[i] <- r[i-1]*ro+x[i]
  s[i]<-  s[i-1]*ro+1
  fg[i] <- s[i-1]*(1-ro)
}
y<-(r+1)/(s+2)


return(list(y,s,fg,ro))

}


# Function that estimates probabilities from a string x
# It returns a list with the estimations, the sample sizes, and the forgotten samples
# It contains a window of active values with fixex size n.


estimate4 <- function(x,n) {
  
  y <- vector(,length(x))  
  s   <- vector(,length(x))  
  ro <-  vector(,length(x))  
  l<-0
  
  k<-1
  
  for(i in 1:length(x)){
    
    if (i-k>=n) {
      l<-1
      k<- k+1
    }
      
    
      
      
    
    
    y[i]<-(sum(x[k:i])+1)/(i-k+1+2)

        s[i] <-  i-k+1
    ro[i] <- l
    
  }
  
  return(list(y,s,ro))
  
}




# Function test with initial simulated data and probability estimations
# It computes the averaged log likelihood of the observations with the estimations of the previous step

test <- function(x,y){
  l <- 0
  n <- length(x)
  for(i in 2:n){
    if(x[i]==1)
    {l <- l + log(y[i-1])} 
    else {l <- l + log(1-y[i-1])}
 
  }

  return(l/(n-1))
  
}


# Function that calls to the appropriate estimation procedure for a set of parameters

sexp <- function(x,param) {
  if(length(param)>0) {
    print(param)
    switch(param[1], 
           '1'=v<-sexp1(x,param), 
           '2'=v<-sexp2(x,param),
           '3'=v<-sexp3(x,param),
           '4'=v<-sexp4(x,param)
           )
    return(v)
  }
  
}

# Function that calls to estimate1 for a set of parameters

sexp1 <- function(x,param) {
  n1 <- as.integer(param[2])
  n2 <- as.integer(param[3])
  h<-sapply(n1:n2, function(y) {z<- estimate1(x,y)
                              plot(z[[1]],type="l")
                             l<- test(x,z[[1]])
                             
                             return(l)
  
  }
)
  met <- rep(1,n2-n1+1)
  arg <- n1:n2
  return(list(h,met,arg))
}

# Function that calls to estimate3 for a set of parameters

sexp3 <- function(x,param) {
  ros <- as.numeric(param[2:length(param)])
  h<-sapply(ros, function(y) {z<- estimate3(x,y)
  plot(z[[1]],type="l")
  l<- test(x,z[[1]])
  
  return(l)
  
  }
  )
  met <- rep(3,length(ros))
  arg <- ros
  return(list(h,met,arg))
}

# Function that calls to estimate2 for a set of parameters


sexp2 <- function(x,param) {
  n1 <- as.integer(param[2])
  n2 <- as.integer(param[3])
  h<-sapply(n1:n2, function(y) {z<- estimate2(x,y)
  plot(z[[1]],type="l")
  l<- test(x,z[[1]])
  return(l)
  }
  )
  met <- rep(2,n2-n1+1)
  arg <- n1:n2
  
  return(list(h,met,arg))
}

# Function that calls to estimate2 for a set of parameters


sexp4 <- function(x,param) {
  n1 <- as.integer(param[2])
  n2 <- as.integer(param[3])
  h<-sapply(n1:n2, function(y) {z<- estimate4(x,y)
  plot(z[[1]],type="l")
  l<- test(x,z[[1]])
  return(l)
  }
  )
  met <- rep(4,n2-n1+1)
  arg <- n1:n2
  
  return(list(h,met,arg))
}

experiment <- function(name){
  
  con <- file(name,"r")
  line <- readLines(con, n = 1)
  param <- strsplit(line,",")
  n <- param[[1]][1]
  nchanges <- param[[1]][2]
  repet <- readLines(con, n = 1)
  met <- readLines(con,n=-1)
close(con)
    mets <- strsplit(met,",") 
   results <- list(vector(),vector(),vector())
  for(i in 1:repet) {
  x <- simulate(n,runif(nchanges))
  for(line in mets) {
    h<- sexp(x,line)
     if(length(h)>0){
        results<- list(c(results[[1]],h[[1]]),  c(results[[2]],h[[2]]), c(results[[3]],h[[3]]))
     }
  }
  }
   names(results) <- c("loglike", "method","param")
   print(results)
   
   mres <- matrix(results[[1]], ncol=as.integer(repet))
   mmethod<- results[[2]][1:nrow(mres)]
   mparam<- results[[3]][1:nrow(mres)]
   
   averages <- rowMeans(mres)
   
  
   print(mmethod)
   print(mparam)
   print(averages)
   return(list(mres,mmethod,mparam))
    
}
experiment("exp1")

sexp3(x,c(3,0.6,0.7,0.8))

x <- simulate(1000,c(0.2,0.5,0.8))
x

t <- estimate3(x,0.9)
t2 <- estimate2(x,30)

plot(t[[1]],type='l')
plot(t2[[2]],type='l')
plot(t2[[3]],type='l')
plot(t2[[4]],type='l')

test(x,t[[1]])
test(x,t2[[1]])

