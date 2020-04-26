library("MASS")
library("abind")
#ģ??????
seed=sample(100000,1)
set.seed(seed)
N=1000
p=10
kk=p+1#??????��
X<-matrix(0,N,kk)
X[,1]=1
S <- toeplitz((p:1)/p)
R <- rWishart(1, p, S)
R=R[,,1]
miu=matrix(0,p,1)
X[,2:(p+1)]=mvrnorm(N,miu,R)
inv=solve(t(X)%*%X)
truepara_1<-rnorm(kk,0,1)
Y=X%*%truepara_1
def=sample(1:N,N/2)##ȱʧ????
Y[def]=NA


INV<-function(x){
  x/sum(x^2)
}
metrics<-function(a,b){
  log((sum((a-b)^2)),10)
}
#####################EM?㷨????############################
oneiter<-function(theta){
  theta=matrix(theta,kk,1)
  #E??
  for(i in def){
    Y[i]=X[i,]%*%theta
  }
  ##M??
  theta=inv%*%t(X)%*%Y
  theta=matrix(theta,1,kk)
  theta
}#һ??????
#####################epsilon_0?????㷨############################
eps_0<-function(theta1,theta2,theta3){
  a=theta3-theta2
  b=theta2-theta1
  theta1-INV(a-b)*sum(b^2)
}##epsilon??????��?㷨???е?һ??????
#####################epsilon_1?????㷨############################
eps_1<-function(theta1,theta2,theta3){
  a=theta3-theta2
  b=theta2-theta1
  theta2+INV(INV(a)-INV(b))
}##epsilon??????��?㷨???е?һ??????
eps<-function(theta1,theta2,theta3){
  theta1+INV(theta3-theta2)
}

###############################һ??ģ???????仯ͼ#####################################
l=300
threshold=-30
thetafirst<-rnorm(kk,0,1)
diff<-matrix(0,l,4)
limit<-matrix(0,l,4)
theta1<-thetafirst#??ʼ????
theta2<-oneiter(theta1)
while(metrics(theta1,theta2)>threshold){
  theta1=theta2
  theta2=oneiter(theta1)
}
truepara=theta2
#########EM(1)
theta<-matrix(0,l,kk)
theta[1,]<-thetafirst#??ʼ????
for(step in 1:(l-1)){
  theta[step+1,]=oneiter(theta[step,])
}
for(step in 2:(l-2)){
  diff[step,1]=metrics(theta[step-1,],theta[step,])
  limit[step,1]=metrics(theta[step,],truepara)
}
#########(2)
theta_1=theta
for(step in 2:(l-1)){
  theta_1[step-1,]=eps_1(theta_1[step-1,],theta_1[step,],theta_1[step+1,])
}
for(step in 2:(l-2)){
  diff[step,2]=metrics(theta_1[step-1,],theta_1[step,])
  limit[step,2]=metrics(theta_1[step,],truepara)
}


###########(3)
cl=ceiling(l/2)
theta_2<-array(rep(0,kk*l*cl),dim=c(l,kk,cl))
theta_2[,,1]<-theta
diff[2,3]=metrics(theta_2[2,,1],theta_2[1,,1])
limit[2,3]=metrics(theta_2[2,,1],truepara)
diff[3,3]=metrics(theta_2[2,,1],theta_2[3,,1])
limit[3,3]=metrics(theta_2[3,,1],truepara)
theta_2[1,,2]=eps_1(theta_2[1,,1],theta_2[2,,1],theta_2[3,,1])
for(step in 4:l){
  diff[step,3]=metrics(theta_2[step,,1],theta_2[step-1,,1])
  limit[step,3]=metrics(theta_2[step,,1],truepara)
  for(k in 2:ceiling(step/2)){
    L=step-2*(k-1)
    theta_2[L,,k]=eps_1(theta_2[L,,k-1],theta_2[L+1,,k-1],theta_2[L+2,,k-1])
    if(L>1){
      diff_1=metrics(theta_2[L,,k],theta_2[L-1,,k])
      limit_1=metrics(theta_2[L,,k],truepara)
      if(is.finite(diff_1)&&diff_1<diff[step,3]){diff[step,3]=diff_1}
      if(is.finite(limit_1)&&limit_1<limit[step,3]){limit[step,3]=limit_1}
    }
  }
}

###########(4)
theta_3<-array(rep(0,kk*l*(l+1)),dim=c(l,kk,l+1))
theta_3[,,2]<-theta
diff[2,4]=metrics(theta_3[2,,2],theta_3[1,,2])
limit[2,4]=metrics(theta_3[2,,2],truepara)
diff[3,4]=metrics(theta_3[2,,2],theta_3[3,,2])
limit[3,4]=metrics(theta_3[3,,2],truepara)
for(j in 3:5){
  for(step in 1:(4-(j-2))){
    theta_3[step,,j]=eps(theta_3[step+1,,j-2],theta_3[step,,j-1],theta_3[step+1,,j-1])
  }
  
}
if(metrics(theta_3[3,,2],theta_3[4,,2])<metrics(theta_3[1,,4],theta_3[2,,4])){
  diff[4,4]=metrics(theta_3[3,,2],theta_3[4,,2])
  limit[4,4]=metrics(theta_3[4,,2],truepara)
}else{
  diff[4,4]=metrics(theta_3[1,,4],theta_3[2,,4])
  limit[4,4]=metrics(theta_3[2,,4],truepara)
}
for(step in 5:l){
  diff[step,4]=metrics(theta_3[step,,2],theta_3[step-1,,2])
  limit[step,4]=metrics(theta_3[step,,2],truepara)
  for(k in 3:(step+1)){
    L=step-k+2
    theta_3[L,,k]=eps(theta_3[L+1,,k-2],theta_3[L,,k-1],theta_3[L+1,,k-1])
    if(k%%2==0&&(L>1)){
      diff_2=metrics(theta_3[L,,k],theta_3[L-1,,k])
      limit_2=metrics(theta_3[L,,k],truepara)
      if(is.finite(diff_2)&&diff_2<diff[step,4]){diff[step,4]=diff_2}
      if(is.finite(limit_2)&&limit_2<limit[step,4]){limit[step,4]=limit_2}
      
    }
  }
}


View(diff[1:(l-4),])
View(limit[1:(l-4),])
plot(diff[2:(l-4),4],type="n",ylab="distance",xlab="k",
     ylim=c(-30,4))
lines(diff[2:(l-4),1],col="red",lty=1,pch=15,lwd=1)
lines(diff[2:(l-4),2],col="green",lty=2,pch=16,lwd=2)
lines(diff[2:(l-4),3],col="blue",lty=3,pch=17,lwd=3)
lines(diff[2:(l-4),4],col="black",lty=4,pch=18,lwd=4)
legend("bottomleft",c("(1)","(2)","(3)","(4)"),
       col=c("red","green","blue","black"),
       text.col=c("red","green","blue","black"),
       pch=c(15,16,17,18),lty=c(1,2,3,4),bty="n",cex=0.6)



plot(limit[2:(l-4),4],type="n",ylab="limit",xlab="k",ylim=c(-30,4))
lines(limit[2:(l-4),1],col="red",lty=1,pch=15,lwd=1)
lines(limit[2:(l-4),2],col="green",lty=2,pch=16,lwd=2)
lines(limit[2:(l-4),3],col="blue",lty=3,pch=17,lwd=3)
lines(limit[2:(l-4),4],col="black",lty=4,pch=18,lwd=4)
legend("bottomleft",c("(1)","(2)","(3)","(4)"),
       col=c("red","green","blue","black"),
       text.col=c("red","green","blue","black"),
       pch=c(15,16,17,18),lty=c(1,2,3,4),bty="n",cex=0.6)



####################################iter_delta2#######################################
start=1
l=20
cl=ceiling(l/2)
theta_4=array(rep(0,kk*l*cl),dim=c(l,kk,cl))
theta_4[,,1]=theta[start:(start+l-1),]
diff4<-matrix(0,l,cl)
for(m in 2:cl){
  for(j in 1:(l-2*(m-1))){
    theta_4[j,,m]=eps_1(theta_4[j,,m-1],theta_4[j+1,,m-1],theta_4[j+2,,m-1])
  }
}
for(j in 1:cl){
  for(i in 1:(l-2*(j-1))){
    diff4[i,j]=metrics(theta_4[i,,j],truepara)
  }
}

View(diff4)

####################################һ????epsilon?㷨#######################################
theta_5=array(rep(0,kk*l*(l+1)),dim=c(l,kk,l+1))
theta_5[,,1]=0
theta_5[,,2]=theta[start:(start+l-1),]
diff5<-matrix(0,l,cl)
for(m in 3:(l+1)){
  for(j in 1:(l-m+2)){
    theta_5[j,,m]=eps(theta_5[j+1,,m-2],theta_5[j,,m-1],theta_5[j+1,,m-1])
  }
}
for(j in 1:cl){
  for(i in 1:(l-2*(j-1))){
    diff5[i,j]=metrics(theta_5[i,,2*j],truepara)
  }
}

View(diff5)








#######################################################################
#####################????ģ??Ч?ʱȽ?##################################
nn=1000
threshold=-10
iternum<-matrix(0,nn,4)##????????
time<-matrix(0,nn,4)##CPUʱ??
for(i in 1:nn){
  thetafirst<-rnorm(kk,0,1)
  #############################em?㷨#################################################
  d=proc.time()
  l=6##??ʼ????????
  theta<-matrix(0,l,kk)
  theta[1,]<-thetafirst
  for(step in 1:(l-1)){
    theta[step+1,]=oneiter(theta[step,])
  }
  diff<-metrics(theta[l,],theta[l-1,])
  while(diff>threshold){
    l=l+1
    theta=rbind(theta,rep(0,kk))
    theta[l,]=oneiter(theta[l-1,])
    diff=metrics(theta[l-1,],theta[l,])
  }
  iternum[i,1]=l
  time[i,1]<-(proc.time()-d)[1]
  
  ##############################eps_1?????㷨#######################################
  d=proc.time()
  l=6
  theta<-matrix(0,l,kk)
  theta[1,]<-thetafirst
  theta[2,]<-oneiter(theta[1,])
  for(step in 2:(l-1)){
    theta[step+1,]=oneiter(theta[step,])
    theta[step-1,]=eps_1(theta[step-1,],theta[step,],theta[step+1,])
  }
  diff=metrics(theta[l-2,],theta[l-3,])
  while(diff>threshold){
    l=l+1
    theta=rbind(theta,rep(0,kk))
    theta[l,]=oneiter(theta[l-1,])
    theta[l-2,]=eps_1(theta[l-2,],theta[l-1,],theta[l,])
    diff=metrics(theta[l-2,],theta[l-3,])
  }
  iternum[i,2]=l
  time[i,2]<-(proc.time()-d)[1]
  
  ############################iteration DELTA2s_1#######################################
  d=proc.time()
  l=3
  cl=ceiling(l/2)
  theta<-array(rep(0,kk*l*cl),dim=c(l,kk,cl))
  theta[1,,1]<-thetafirst
  theta[2,,1]<-oneiter(theta[1,,1])
  theta[3,,1]<-oneiter(theta[2,,1])
  theta[1,,2]=eps_1(theta[1,,1],theta[2,,1],theta[3,,1])
  diff=metrics(theta[3,,1],theta[2,,1])
  while(diff>threshold){
    l=l+1
    cl=ceiling((l-1)/2)
    theta=abind(theta,array(rep(0,cl*kk),dim=c(1,kk,cl)),along=1)
    if(l%%2==1){
      theta=abind(theta,array(rep(0,l*kk),dim=c(l,kk,1)),along=3)
    }
    theta[l,,1]=oneiter(theta[l-1,,1])
    for(k in 2:cl){
      L=l-2*(k-1)
      theta[L,,k]=eps_1(theta[L,,k-1],theta[L+1,,k-1],theta[L+2,,k-1])
      if(L>1){
        diff_1=metrics(theta[L,,k],theta[L-1,,k])
        if(is.finite(diff_1)&&diff_1<diff){diff=diff_1;result_1=theta[L,,k]}
      }
      
    }
  }
  
  iternum[i,3]=l
  time[i,3]<-(proc.time()-d)[1]
  ############################general epsilon#######################################
  d=proc.time()
  l=4
  theta<-array(rep(0,kk*l*(l+1)),dim=c(l,kk,l+1))
  theta[,,1]=0
  theta[1,,2]<-thetafirst
  for(step in 1:(l-1)){
    theta[step+1,,2]=oneiter(theta[step,,2])
  }
  for(j in (l-1):(l+1)){
    for(step in 1:(l-(j-2))){
      theta[step,,j]=eps(theta[step+1,,j-2],theta[step,,j-1],theta[step+1,,j-1])
    }
  }
  if(metrics(theta[1,,4],theta[2,,4])<metrics(theta[3,,2],theta[4,,2])){
    diff=metrics(theta[1,,4],theta[2,,4])
  }else{
    diff=metrics(theta[3,,2],theta[4,,2])
  }
  while(diff>threshold){
    l=l+1
    theta=abind(theta,array(rep(0,kk*(l-1)),dim=c(l-1,kk,1)),along=3)
    theta=abind(theta,array(rep(0,(l+1)*kk),dim=c(1,kk,l+1)),along=1)
    theta[l,,1]=0
    for(k in 2:(l+1)){
      L=l-(k-2)
      if(k==2){
        theta[l,,2]=oneiter(theta[l-1,,2])
        diff=metrics(theta[l,,2],theta[l-1,,2])
        result_2=theta[l,,2]
      }else{
        theta[L,,k]=eps(theta[L+1,,k-2],theta[L,,k-1],theta[L+1,,k-1])
        if(k%%2==0&&(L>1)){
          diff_2=metrics(theta[L,,k],theta[L-1,,k])
          if(is.finite(diff_2)&&diff_2<diff){diff=diff_2;result_2=theta[L,,k]}
        }
      }
    }
  }
  iternum[i,4]=l
  time[i,4]<-(proc.time()-d)[1]
}

par(mfrow = c(2, 4))
hist(iternum[,1],freq=TRUE,xlab="k",ylab="frequency",main="(1)",col="red")
hist(iternum[,2],freq=TRUE,xlab="k",ylab="frequency",main="(2)",col="green")
hist(iternum[,3],freq=TRUE,xlab="k",ylab="frequency",main="(3)",col="blue")
hist(iternum[,4],freq=TRUE,xlab="k",ylab="frequency",main="(4)",col="black")
hist(time[,1],freq=TRUE,xlab="k",ylab="time",main="(1)",col="red")
hist(time[,2],freq=TRUE,xlab="k",ylab="time",main="(2)",col="green")
hist(time[,3],freq=TRUE,xlab="k",ylab="time",main="(3)",col="blue")
hist(time[,4],freq=TRUE,xlab="k",ylab="time",main="(4)",col="black")


summary(iternum)
summary(time)








