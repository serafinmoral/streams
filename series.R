
x1<- rbinom(1000,1,0.2)
x2<- rbinom(1000,1,0.5)
x3<- rbinom(1000,1,0.8)
x <- c(x1,x2,x3)
y<-x
for(i in 1:3000){
    y[i] <- (sum(x[1:i])+1)/(i+2.0)}
y
plot(y,type='l')

s <- vector(,3000)
r <- vector(,3000)
ro  <- vector(,3000)
ros  <- vector(,3000)
ros <- 0.999

y2 <- vector(,3000)
r[1] = x[1]
s[1] = 1
for(i in 2:3000){
  r[i] <- r[i-1]*ros+x[i]
  s[i]<-  s[i-1]*ros+1
  }
y2<-(r+1)/(s+2)
plot(y2,type='l')
y2
ro=1.000
s



y3 <- vector(,3000)
r[1] = x[1]
s[1] = 1
y3[1] <- (r[1]+1)/(s[1]+2.0)
ro[1]<- 1
for(i in 2:3000){
  j<-max(c(1,i-29))
  pvalue <- binom.test(sum(x[j:i]),i-j+1,y3[j])$p.value
  if(pvalue>0.05) {ro[i]<-1}
  else {ro[i] <- pvalue^{1/200}}
  
  r[i] <- r[i-1]*ro[i]+x[i]
  s[i]<-  s[i-1]*ro[i]+1
  y3[i]<- (r[i]+1)/(s[i]+2.0)
}
y3<-(r+1)/(s+2)
plot(y3,type='l')
plot(ro,type='l')
plot(s,type='l')
plot(ros,type='l')
s
ro[1]<- 1
r[1] = x[1]
s[1] = 1

r[1] = x[1]
s[1] = 1
y3[1] <- (r[1]+1)/(s[1]+2.0)
ro[1]<- 1
ros[1]<- 1
for(i in 2:3000){
  j<-max(c(1,i-29))
  pvalue <- binom.test(sum(x[j:i]),i-j+1,y3[j])$p.value
  if(pvalue>0.05) {ro[i]<-1}
  else {ro[i] <- pvalue^{1/200}}
   ros[i]=sum(ro[j:i])/(i-j+1.0) 
  r[i] <- r[i-1]*ros[i]+x[i]
  s[i]<-  s[i-1]*ros[i]+1
  y3[i]<- (r[i]+1)/(s[i]+2.0)
}

for(i in 2:3000){
  if (i<60) {
      ro[i]= 1.0
  }
  else{
  change <-2*lgamma(2) - 2*lgamma(32) - 2*lgamma(1) +
           lgamma(sum(x[(i-30+1):i]) +1 ) + lgamma(sum(x[(i-60+1):(i-30)]) +1) + 
    lgamma(30-sum(x[(i-30+1):i]) +1 ) + lgamma(30- sum(x[(i-60+1):(i-30)]) +1) 
  nochange <-  lgamma(4) - lgamma(64) - 2*lgamma(2) + lgamma(sum(x[(i-60+1):i]) +2) +lgamma(60- sum(x[(i-60+1):i]) +2)
  relative <- nochange-change
  odd <- exp(relative)
  if (odd>0.5) {ro[i]<- 1.0}
  else {ro[i]<- odd^{1/20}}}
  r[i] <- r[i-1]*ro[i]+x[i]
  s[i]<-  s[i-1]*ro[i]+1
  y3[i]<- (r[i]+1)/(s[i]+2.0)
}

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


k<-1

for(i in 2:3000){
  
  l<-0
  if (i-k>=100) {
    test<- 1
    
    while(test==1) {
      j = k+(i-k)%/%2  
    
      change <-(lgamma(2) - lgamma(i-j+1+2)     + 
                  lgamma(sum(x[j:i]) +1 ) -  lgamma(1)   +
                  lgamma((i-j+1-sum(x[j:i]) +1) ) -  lgamma(1)
                + lgamma(2) - lgamma(j-k+2) +
                   lgamma(sum(x[k:(j-1)]) +1 ) -  lgamma(1)   +
                  lgamma((j-k-sum(x[k:(j-1)]) +1) ) -  lgamma(1)
      )
      
      nochange <-  (lgamma(4) - lgamma(i-k+4+1) 
                    - 2*lgamma(2) + 
                      lgamma(sum(x[k:i]) +2) +lgamma(i-k+1 - sum(x[k:i]) +2)
      )
      
      relative <- nochange-change
      
      print(change)
      print(nochange)
      
      odd <- exp(relative) 
      test<-0
      if (odd <0.1) {k<-j}
    }
  }
  
  
  y3[i]<-(sum(x[k:i])+1)/(i-k+1+2)
  s[i] <-  i-k+1
  ro[i] = l
  
}

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

