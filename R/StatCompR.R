 fr=function(x){
    z=rnorm(100)
    d=sort(abs(z-x))
    y=numeric(100)
    for(k in 1:100){
      y[k]=k/2/100/d[k]
    }
    dis=abs(y-dnorm(x))
    re=which.min(dis)
    return(re)
  }
