



ro[1]<- 1
r[1] = x[1]
s[1] = 1




plot(y3,type='l')
plot(ro,type='l')
plot(s,type='l')

k<-1

for(i in 2:3000){
   
    l<-0
  if (i-k>=59) {
    test<- 1
    k
     while(test==1) {
     
       change <-(lgamma(2) - lgamma(32)     +  
                   lgamma(sum(x[(i-30+1):i]) +1 )  -  lgamma(1)   
                 + lgamma((30-sum(x[(i-30+1):i]) +1) ) -  lgamma(1) + 
                   lgamma(2) - lgamma(i-30-k+1+2)
     + lgamma(sum(x[k:(i-30)]) +1) -  lgamma(1)
    +   lgamma(i-30+1-k-sum(x[k:(i-30)]) +1 ) -lgamma(1)
       )
     
     nochange <-  (lgamma(4) - lgamma(i-k+4+1) 
     - 2*lgamma(2) + 
       lgamma(sum(x[k:i]) +2) +lgamma(i-k+1 - sum(x[k:i]) +2)
     )
     
    relative <- nochange-change
  
    
    odd <- exp(relative) 
 test <- 0
if (odd <0.1) {k<- k + (i-30-k)*(1-odd^4)}
     }
  }
    
  
    y3[i]<-(sum(x[k:i])+1)/(i-k+1+2)
    s[i] <-  i-k+1
    ro[i] = l
    
  }

plot(y3,type='l')
plot(ro,type='l')
plot(s,type='l')



plot(y3,type='l')
plot(ro,type='l')
plot(s,type='l')

r2 <- x[1];
s2 <- 1;


for (i in 1:3000) {
  x[i] =  rbinom(1,1, (i/2999-1/2999))
}
y3[2950:3000]




y3 <- vector(,3000)
r2 <- vector(,3000)
s2  <- vector(,3000)
r[1] = x[1]
s[1] = 1
y3[1] <- (r[1]+1)/(s[1]+2.0)
ro[1]<- 1
r2[1] =x[1]
s2[1]=1

for(i in 2:3000){
  if(i<31) {ro[i]<-1
      r[i] <- r[i-1]+x[i]
      s[i] <-s[i-1] + 1
  }
  else {
    j <- i-29
    pvalue <- binom.test(sum(x[j:i]),30,y3[j])$p.value 
    if(pvalue>0.01) {ro[i]<-1}
    else {ro[i] <- pvalue^{1/100}}
    if(i>31) {
    r2[i-30] = r2[i-31]*ro[i] + x[i-30]
    s2[i-30] = s2[i-31]*ro[i] + 1}
    r[i] = r2[i-30] + sum(x[j:i])
    s[i] = s2[i-30] + 30
  }

  y3[i]<- (r[i]+1)/(s[i]+2)
}

plot(y3,type='l')
plot(ro,type='l')
plot(s,type='l')
plot(ros,type='l')

