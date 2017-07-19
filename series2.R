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

#simulate a series of n binomial values with a binomial probability
# which is linearly changing

simulate2 <- function(n) {
  x <- c()
  
  l1 <- runif(1)
  l2 <- runif(1)
  
  
  
   for(i in 1:n)  {
    x1<- rbinom(1,1, l1+(l2-l1)*i/n)
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
# alpha is the threshold for the odds

estimate1 <- function(x,n,alpha) {
  
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
      
      if (odd < alpha) {k<-k+n
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
# m is a value such that the p-value is squared to m

estimate2 <- function(x,n1,m){

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
  else {ro[i] <- pvalue^{1/m}}
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


# Function that estimates the probabilities with a  forgetting factor
# which is computed by a Bayesian BDEu score of the probability of change
# in the last n cases (dividing these cases in two equal parts)


estimate5 <- function(x,n) {
  
  y <- vector(,length(x))  
  s   <- vector(,length(x))  
  ro <- rep(1,length(x))
  fg <- vector(,length(x)) 
  r <- vector(,length(x))  
  
  r[1] = x[1]
  s[1] = 1
  fg[1]=0
  ro[1] = 1
  y[1]<- (r[1]+1)/(s[1]+2.0)
  
for(i in 2:length(x)){
  if (i<n) {
    ro[i]= 1.0
  }
  else{
    n1 <- n %/% 2
    n2 <- n-n1
     odd <- computeoddsc(n1,n2,sum(x[(i-n+1):(i-n2)]), sum(x[(i-n2+1):(i)]),4  )
    if (odd>alpha) {ro[i]<- 1.0}
    else {ro[i]<- odd^{1/20}}}
  r[i] <- r[i-1]*ro[i]+x[i]
  s[i]<-  s[i-1]*ro[i]+1
  fg[i] = s[i-1] * (1-ro[i])
  y[i]<- (r[i]+1)/(s[i]+2.0)
}

  
  
  return(list(y,s,fg,ro))
  
}


# Function that estimates probabilities from a string x
# It returns a list with the estimations, the sample sizes, and the forgotten samples
# It contains a window of active values. 
# It is similar to estimate1, but now the window of active 
# values is divided by two, if the total size is greater
# than n


estimate6 <- function(x,n,alpha) {
  
  y <- vector(,length(x))  
  s   <- vector(,length(x))  
  ro <-  vector(,length(x))  
  l<-0
  
  k<-1
  
  for(i in 1:length(x)){
    
    
    
    if (i-k>=2*n-1) {
      
      l<- 0
      
        j1 = k+ (i-k)%/%2
        j2 = j1+1
        
    
        x1 = sum(x[k:j1])
        x2 = sum(x[j2:i])
        
        odd <- computeoddsc(j1-k+1,i-j2+1,x1,x2,4)
        
        if (odd <alpha) { l<- j1-k+1
        k<-j2
        
        
        }
      
      }
      
    
    
    y[i]<-(sum(x[k:i])+1)/(i-k+1+2)
    s[i] <-  i-k+1
    ro[i] <- l
    
  }
  
  return(list(y,s,ro))
  
}


# Function that estimates probabilities from a string x
# It returns a list with the estimations, the sample sizes, and the forgotten samples
# It contains a window of active values. 
# It is similar to estimate1 and estimate5, but now the window of active values is divided in two
# parts if its length is greater than n: one with the last n values and other with the rest of values


estimate7 <- function(x,n,alpha) {
  
  y <- vector(,length(x))  
  s   <- vector(,length(x))  
  ro <-  vector(,length(x))  
  l<-0
  
  k<-1
  
  for(i in 1:length(x)){
    
    
    
    if (i-k>=n) {
      
      l<- 0
      
      j1 = i-n
      j2 = i-n+1
      
      
      x1 = sum(x[k:j1])
      x2 = sum(x[j2:i])
      
      odd <- computeoddsc(i-k-n+1,n,x1,x2,4)
      
      if (odd <alpha) {
        l<- j1-k+1
        k<-j2
     
      
      }
      
    }
    
    
    
    y[i]<-(sum(x[k:i])+1)/(i-k+1+2)
    s[i] <-  i-k+1
    ro[i] <- l
    
  }
  
  return(list(y,s,ro))
  
}


# Function that estimates probabilities from a string x
# It returns a list with the estimations, the sample sizes, and the forgotten samples
# It contains a window of active values. 
# It is similar to estimate7, but now it carries out a chisquared test


estimate8 <- function(x,n,alpha) {
  
  require(MASS)
  y <- vector(,length(x))  
  s   <- vector(,length(x))  
  ro <-  vector(,length(x))  
  l<-0
  
  k<-1
  
  for(i in 1:length(x)){
    
    
    
    if (i-k>=n) {
      
      l<- 0
      
      j1 = i-n
      j2 = i-n+1
      
      
      x1 = sum(x[k:j1])
      x2 = sum(x[j2:i])
      x1n = j1-k+1 - x1
      x2n <- i-j2+1 - x2
      
      total <- i-k+1
      
      e11 <- (x1+x1n)*(x1+x2)/total
      e12 <-  (x2+x2n)*(x1+x2)/total
      e21 <-  (x1+x1n)*(x1n+x2n)/total
      e22 <-  (x2+x2n)*(x1n+x2n)/total
      
      tab <- array(c(x1,x1n,x2,x2n), dim = c(2,2))
      
      if((e11>=5) && (e12>=5) && (e21>=5) && (e22>=5)) {
        p<-   chisq.test(tab)$p.value
      }
      else 
      {
        p <- fisher.test(tab)$p.value
      }
      
      
      if (p <alpha) {
        l<- j1-k+1
        k<-j2
        
        
      }
      
    }
    
    
    
    y[i]<-(sum(x[k:i])+1)/(i-k+1+2)
    s[i] <-  i-k+1
    ro[i] <- l
    
  }
  
  return(list(y,s,ro))
  
}


# Function that estimates probabilities from a string x
# It returns a list with the estimations, the sample sizes, and the forgotten samples
# forg is a vector of rhos. It considers all the rhos and makes an average of the estimations
# with the different rhos weighted by the likelihood of the data given each rho


estimate9<- function(x,forg,l) {
  count <- vector("double",l)
  sample <- vector("double",l)
  logd <-  vector("double",l)
  
  
  y <- vector("double",length(x))  
  s   <- vector("double",length(x))  
  ro <-  vector("double",length(x))  
  
  count = rep(x[1],l)
  sample = rep(1,l)
  
  
  
  
  y[1] <- (x[1]+1)/3
  s[1] = 1
  ro[1] = 1
  
  for(i in 2:length(x)){
      for(j in 1:l) {
  
        count[j] <- x[i] + forg[j] * count[j]
        sample[j] <- 1 + forg[j] * sample[j]
  #     print(count[j])
        if(y[1] ==1) {logd[j] <- logd[j] + log((count[j]+1)/(sample[j]+2))}
        else {logd[j] <-logd[j]  +log( (sample[j]-count[j]+1)/(sample[j]+2))}
      }
    logd = logd - max(logd)
    y[i] <-0
    s[i] <-0
    sumt <- sum(exp(logd))
    print(sumt)
    for(j in 1:l) {
      y[i] <- y[i] + exp(logd[j]) * (count[j]+1)/(sample[j]+2)
      s[i] <- s[i] + exp(logd[j]) *(sample[j])   
      ro[i] =  ro[i] + exp(logd[j]) *forg[j] 
    }
    y[i]<- y[i]/sumt
    print(y[i])
    s[i] <- s[i]/sumt
    ro[i] <- ro[i]/sumt
    
  }
  
  
  return(list(y,s,ro)) 
  
}


# Function that estimates probabilities from a string x
# It returns a list with the estimations, the sample sizes, and the forgotten samples
# forg is a vector of rhos. It considers all the rhos and selects 
# the rho with maximum likelihood in each case.


estimate10<- function(x,forg,l) {
  count <- vector("double",l)
  sample <- vector("double",l)
  logd <-  vector("double",l)
  
  
  y <- vector("double",length(x))  
  s   <- vector("double",length(x))  
  ro <-  vector("double",length(x))  
  
  count = rep(x[1],l)
  sample = rep(1,l)
  
  
  
  
  y[1] <- (x[1]+1)/3
  s[1] = 1
  ro[1] = 1
  
  for(i in 2:length(x)){
    for(j in 1:l) {
      
      count[j] <- x[i] + forg[j] * count[j]
      sample[j] <- 1 + forg[j] * sample[j]
      #     print(count[j])
      if(y[1] ==1) {logd[j] <- logd[j] + log((count[j]+1)/(sample[j]+2))}
      else {logd[j] <-logd[j]  +log( (sample[j]-count[j]+1)/(sample[j]+2))}
    }
    logd = logd - max(logd)
    j <- which.max(logd)
    
    
  
      y[i] <- (count[j]+1)/(sample[j]+2)
      s[i] <- (sample[j])   
      ro[i] =  forg[j] 
    
      
    
    
  }
  
  
  return(list(y,s,ro)) 
  
}

# Function that estimates probabilities according the Bifet, Gavalda, 2007 procedure


estimate11<- function(x,delta) {
  
  
  y <- vector("double",length(x))  
  s   <- vector("double",length(x))  
  ro <-  vector("double",length(x))  
  
  
  y[1] <- (x[1]+1)/3
  s[1] = 1
  ro[1] = 1
  k<- 1
  
  for(i in 2:length(x)){
    fin<- FALSE
    while(!fin){
      fin <- TRUE
    for(j in k:(i-1)) {
       n1 <- sum(x[k:j])
       n2 <- sum(x[(j+1):i])
       l1 <- j-k+1
       l2 <- i-j
       
       m <- 1/(1/l1 + 1/l2)
       deltap <- delta/log(i-k+1)
       cut<- sqrt(1/(2*m) *( (n1+n2)/(i-k+1))*
                    ((i-k+1 -n1-n2)/(i-k+1)) * 
                    log(2/deltap)) + 2/(3*m) * 
         log(2/deltap)
       if (abs(n1/l1-n2/l2)> cut) {
         fin <- FALSE
         k <- k+1
         break
       }
    }
  }
    y[i] <- (sum(x[k:i])+1)/(i-k+3)
    s[i]<- i-k+1
    ro[i] <- (s[i]-1)/s[i-1]
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
           '4'=v<-sexp4(x,param),
           '5'=v<-sexp5(x,param),
           '6'=v<-sexp6(x,param),
           '7'=v<-sexp7(x,param),
           '8'=v<-sexp8(x,param),
           '9'=v<-sexp9(x,param),
           '10'=v<-sexp10(x,param)
           
           )
    return(v)
  }
  
}

# Function that calls to estimate1 for a set of parameters

sexp1 <- function(x,param) {
  n1 <- as.integer(param[2])
  n2 <- as.integer(param[3])
  alpha <- as.numeric(param[4])
  h<-sapply(n1:n2, function(y) {
                          z<- estimate1(x,y,alpha)
                              plot(z[[1]],type="l")
                             l<- test(x,z[[1]])
                             
                             return(l)
  
  }
)
  met <- rep(1,n2-n1+1)
  arg1 <- n1:n2
  arg2 <- rep(alpha,n2-n1+1)
  
  return(list(h,met,arg1,arg2))
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
  h<-sapply(n1:n2, function(y) {z<- estimate2(x,y,100)
  plot(z[[1]],type="l")
  l<- test(x,z[[1]])
  return(l)
  }
  )
  met <- rep(2,n2-n1+1)
  arg <- n1:n2
  
  return(list(h,met,arg))
}

# Function that calls to estimate4 for a set of parameters


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



# Function that calls to estimate5 for a set of parameters


sexp5 <- function(x,param) {
  n1 <- as.integer(param[2])
  n2 <- as.integer(param[3])
  h<-sapply(n1:n2, function(y) {z<- estimate5(x,y,0.5)
  plot(z[[1]],type="l")
  l<- test(x,z[[1]])
  return(l)
  }
  )
  met <- rep(5,n2-n1+1)
  arg <- n1:n2
  
  return(list(h,met,arg))
}


# Function that calls to estimate6 for a set of parameters


sexp6 <- function(x,param) {
  n1 <- as.integer(param[2])
  n2 <- as.integer(param[3])
  h<-sapply(n1:n2, function(y) {z<- estimate6(x,y,0.5)
  plot(z[[1]],type="l")
  l<- test(x,z[[1]])
  return(l)
  }
  )
  met <- rep(6,n2-n1+1)
  arg <- n1:n2
  
  return(list(h,met,arg))
}



# Function that calls to estimate7 for a set of parameters


sexp7 <- function(x,param) {
  n1 <- as.integer(param[2])
  n2 <- as.integer(param[3])
  h<-sapply(n1:n2, function(y) {z<- estimate7(x,y,0.1)
  plot(z[[1]],type="l")
  l<- test(x,z[[1]])
  return(l)
  }
  )
  met <- rep(7,n2-n1+1)
  arg <- n1:n2
  
  return(list(h,met,arg))
}



# Function that calls to estimate8 for a set of parameters


sexp8 <- function(x,param) {
  n1 <- as.integer(param[2])
  n2 <- as.integer(param[3])
  h<-sapply(n1:n2, function(y) {z<- estimate8(x,y,0.005)
  plot(z[[1]],type="l")
  l<- test(x,z[[1]])
  return(l)
  }
  )
  met <- rep(8,n2-n1+1)
  arg <- n1:n2
  
  return(list(h,met,arg))
}


# Function that calls to estimate9 for a set of parameters


sexp9 <- function(x,param) {
   l <- length(param)-1
  
   
   forg <- as.double(param[2:(l+1)])
   
   z<- estimate9(x,forg,l)
   plot(z[[1]],type="l")
   l<- test(x,z[[1]])
   par <- list(forg)
   return(list(l,9,par))
}



# Function that calls to estimate10 for a set of parameters


sexp10 <- function(x,param) {
  l <- length(param)-1
  
  
  forg <- as.double(param[2:(l+1)])
  
  z<- estimate10(x,forg,l)
  plot(z[[1]],type="l")
  l<- test(x,z[[1]])
  par <- list(forg)
  return(list(l,10,par))
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



experiment2 <- function(name){
  
  con <- file(name,"r")
  line <- readLines(con, n = 1)
  n <- as.integer(line)
  repet <- readLines(con, n = 1)
  met <- readLines(con,n=-1)
  close(con)
  mets <- strsplit(met,",") 
  results <- list(vector(),vector(),vector())
  for(i in 1:repet) {
    x <- simulate2(n)
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
experiment2("exp2")

sexp5(x,c(3,100,101,102))

x <- simulate(1000,c(0.2,0.5,0.8))
x

t <- estimate8(x,50)
t2 <- estimate3(x,0.99)

t3 <- estimate4(x,200)

test(x,t[[1]])

test(x,t2[[1]])

test(x,t3[[1]])


plot(t[[1]],type='l')
plot(t2[[1]],type='l')
plot(t3[[1]],type='l')
plot(t2[[3]],type='l')

test(x,t[[1]])
test(x,t2[[1]])

