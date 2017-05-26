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
  ro[i] = l
  
}

return(list(y,s,ro))

}

test <- function(x,y){
  l <- 0
  n <- length(x)
  for(i in 2:n){
    l <- l + log(y[i-1])^x[i] + log(1-y[i-1])^(1-x[i])
  }
  return(l/n)
  
}


x <- simulate(1000,c(0.2,0.5,0.8))
t <- estimate1(x,200)

plot(t[[1]],type='l')
plot(t[[2]],type='l')
plot(t[[3]],type='l')

t[1]
