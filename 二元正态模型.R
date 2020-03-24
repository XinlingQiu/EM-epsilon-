#####################################################
#模拟样本
library('abind')
seed=sample(100000,1)
set.seed(seed)
x1<-c(8,11,16,18,25,9,13,NaN,NaN,NaN);
x2<-c(10,14,16,15,NaN,NaN,NaN,15,20,4);
N=10;
kk=5
k=1;
dex=0;
dey=0;
for(i in 1:N){
  if(is.nan(x1[i])){dex[k]=i;k=k+1;}
}
k=1;
for(i in 1:N){
  if(is.nan(x2[i])){dey[k]=i;k=k+1;}
}


INV<-function(x){
  x/sum(x^2)
}
metrics<-function(a,b){
  log((sum((a-b)^2)),10)
}
#####################EM算法序列############################
oneiter<-function(theta){
  miu=theta[1:2]
  sig=theta[3:5]
  #E步
  for(i in dex){
    x1[i]=miu[1]+sig[3]*(x2[i]-miu[2])/sig[2];
  }
  for(i in dey){
    x2[i]=miu[2]+sig[3]*(x1[i]-miu[1])/sig[1];
  }
  #M步
  m1=sum(x1);
  m2=sum(x2);
  miu[1]=m1/N;
  miu[2]=m2/N;
  sig[1]=sum((x1-miu[1])^2)/N;
  sig[2]=sum((x2-miu[2])^2)/N;
  sig[3]=sum(t(x1-miu[1])%*%(x2-miu[2]))/N;
  theta=c(miu,sig)
  theta
}#一步迭代
#####################epsilon_0加速算法############################
eps_0<-function(theta1,theta2,theta3){
  a=theta3-theta2
  b=theta2-theta1
  theta1-INV(a-b)*sum(b^2)
}##epsilon加速收敛算法序列的一步迭代
#####################epsilon_1加速算法############################
eps_1<-function(theta1,theta2,theta3){
  a=theta3-theta2
  b=theta2-theta1
  theta2+INV(INV(a)-INV(b))
}##epsilon加速收敛算法序列的一步迭代
#####################epsilon_2加速算法############################
eps_2<-function(theta1,theta2,theta3){
  a=theta3-theta2
  b=theta2-theta1
  theta3-INV(a-b)*sum(a^2)
}##epsilon加速收敛算法序列的一步迭代

#####################iteration delta2_0算法############################
rdel2_0<-function(theta){
  l=nrow(theta)
  while(l-2>0){
    for(j in 1:(l-2)){
      theta[j,]=eps_0(theta[j,],theta[j+1,],theta[j+2,])
    }
    l=l-2
  }
  theta
}
#####################iteration delta2_1算法############################
rdel2_1<-function(theta){
  l=nrow(theta)
  while(l-2>0){
    for(j in 1:(l-2)){
      theta[j,]=eps_1(theta[j,],theta[j+1,],theta[j+2,])
    }
    l=l-2
  }
  theta
}
#####################iteration delta2_2算法############################
rdel2_2<-function(theta){
  l=nrow(theta)
  while(l-2>0){
    for(j in 1:(l-2)){
      theta[j,]=eps_2(theta[j,],theta[j+1,],theta[j+2,])
    }
    l=l-2
  }
  theta
}
#####################一般性epsilon算法完全迭代############################
geps<-function(theta){
  l=nrow(theta)
  t<-matrix(0,l,kk)
  l=l-1
  while(l>0){
    for(k in 1:l){
      t[k,]=t[k+1,]+INV(theta[k+1,]-theta[k,])
    }
    l=l-1
    if(l>0){
      for(k in 1:l){
        theta[k,]=theta[k+1,]+INV(t[k+1,]-t[k,])
      }
      l=l-1
    }
  }
  theta
}##一般性epsilon加速收敛算法序列的完全迭代


###############################一次模拟误差变化图#####################################
l=200
threshold=-30
miufirst<-runif(2)
sigfirst<-runif(3)
thetafirst<-c(miufirst,sigfirst)

diff<-matrix(0,l,4)
limit<-matrix(0,l,4)
theta1<-thetafirst#初始参数
theta2<-oneiter(theta1)
while(metrics(theta1,theta2)>threshold){
  theta1=theta2
  theta2=oneiter(theta1)
}
truepara=theta2
#truepara=c(13.673,13.959 ,53.017 ,22.061 ,32.910)
#########EM
theta<-matrix(0,l,kk)
theta[1,]<-thetafirst#初始参数
for(step in 1:(l-1)){
  theta[step+1,]=oneiter(theta[step,])
}
for(step in 2:(l-2)){
  diff[step,1]=metrics(theta[step-1,],theta[step,])
  limit[step,1]=metrics(theta[step,],truepara)
}


#########epsilon_1加速
theta_1=theta
for(step in 2:(l-1)){
  theta_1[step-1,]=eps_1(theta_1[step-1,],theta_1[step,],theta_1[step+1,])
}
for(step in 2:(l-2)){
  diff[step,2]=metrics(theta_1[step-1,],theta_1[step,])
  limit[step,2]=metrics(theta_1[step,],truepara)
}


###########rdel2_1
theta_2=theta_1
for(step in 2:(l-3)){
  theta_2[step-1,]=eps_1(theta_2[step-1,],theta_2[step,],theta_2[step+1,])
}
for(step in 2:(l-4)){
  diff[step,3]=metrics(theta_2[step-1,],theta_2[step,])
  limit[step,3]=metrics(theta_2[step,],truepara)
}
###########一般性epsilon
theta_6<-array(rep(0,kk*l*(l+1)),dim=c(l,kk,l+1))
theta_6[,,2]<-theta
diff[2,4]=metrics(theta_6[2,,2],theta_6[1,,2])
limit[2,4]=metrics(theta_6[2,,2],truepara)
diff[3,4]=metrics(theta_6[2,,2],theta_6[3,,2])
limit[3,4]=metrics(theta_6[3,,2],truepara)
for(j in 3:5){
  for(step in 1:(4-(j-2))){
    theta_6[step,,j]=theta_6[step+1,,j-2]+INV(theta_6[step+1,,j-1]-theta_6[step,,j-1])
  }
}
if(metrics(theta_6[3,,2],theta_6[4,,2])<metrics(theta_6[1,,4],theta_6[2,,4])){
  diff[4,4]=metrics(theta_6[3,,2],theta_6[4,,2])
  limit[4,4]=metrics(theta_6[4,,2],truepara)
}else{
  diff[4,4]=metrics(theta_6[1,,4],theta_6[2,,4])
  limit[4,4]=metrics(theta_6[2,,4],truepara)
}
for(step in 5:l){
  diff[step,4]=metrics(theta_6[step,,2],theta_6[step-1,,2])
  limit[step,4]=metrics(theta_6[step,,2],truepara)
  for(k in 3:(step+1)){
    L=step-k+2
    theta_6[L,,k]=theta_6[L+1,,k-2]+INV(theta_6[L+1,,k-1]-theta_6[L,,k-1])
    if(k%%2==0){
      if(step%%2==0||(step%%2==1&&k<step)){
        diff_1=metrics(theta_6[L,,k],theta_6[L-1,,k])
        limit_1=metrics(theta_6[L,,k],truepara)
        if(diff_1<diff[step,4]){diff[step,4]=diff_1}
        if(limit_1<limit[step,4]){limit[step,4]=limit_1}
      }
    }
  }
}


View(diff[1:(l-4),])
View(limit[1:(l-4),])
plot(diff[2:(l-4),4],type="n",ylab="收敛误差",xlab="迭代次数",
     ylim=c(-30,4))
lines(diff[2:(l-4),1],col="red",lty=1,pch=15,lwd=1)
lines(diff[2:(l-4),2],col="green",lty=2,pch=16,lwd=2)
lines(diff[2:(l-4),3],col="blue",lty=3,pch=17,lwd=3)
lines(diff[2:(l-4),4],col="black",lty=4,pch=18,lwd=4)
legend("bottomleft",c("em","eps","迭代delta2",
                      "一般性epsilon"),
       col=c("red","green","blue","black"),
       text.col=c("red","green","blue","black"),
       pch=c(15,16,17,18),lty=c(1,2,3,4),bty="n",cex=0.6)



plot(limit[2:(l-4),4],type="n",ylab="与真值的误差",xlab="迭代次数",ylim=c(-30,4))
lines(limit[2:(l-4),1],col="red",lty=1,pch=15,lwd=1)
lines(limit[2:(l-4),2],col="green",lty=2,pch=16,lwd=2)
lines(limit[2:(l-4),3],col="blue",lty=3,pch=17,lwd=3)
lines(limit[2:(l-4),4],col="black",lty=4,pch=18,lwd=4)
legend("bottomleft",c("em","eps","迭代delta2",
                      "一般性epsilon"),
       col=c("red","green","blue","black"),
       text.col=c("red","green","blue","black"),
       pch=c(15,16,17,18),lty=c(1,2,3,4),bty="n",cex=0.6)



####################################一般性epsilon算法#######################################
l=100
cl=ceiling(l/2)
#############一般性epsilon算法完全迭代
theta_3=array(rep(0,kk*l*(l+1)),dim=c(l,kk,l+1))
theta_3[,,1]=0
theta_3[,,2]=theta
diff3<-matrix(0,l,cl)
for(m in 3:(l+1)){
  for(j in 1:(l-m+2)){
    theta_3[j,,m]=theta_3[j+1,,m-2]+INV(theta_3[j+1,,m-1]-theta_3[j,,m-1])
  }
}
for(j in 1:cl){
  for(i in 1:(l-2*(j-1))){
    diff3[i,j]=metrics(theta_3[i,,2*j],truepara)
  }
}

View(diff3)








#######################################################################
#####################多次模拟效率比较##################################
nn=1000
threshold=-10
iternum<-matrix(0,nn,4)##迭代次数
time<-matrix(0,nn,4)##CPU时间
theta<-array(rep(0,kk*l*cl),dim=c(l,kk,cl))
for(i in 1:nn){
  miufirst<-runif(2)
  sigfirst<-runif(3)
  thetafirst<-c(miufirst,sigfirst)
  
  #############################em算法#################################################
  d=proc.time()
  l=6##初始迭代次数
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

  ##############################eps_1加速算法#######################################
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
  l=6
  theta<-array(rep(0,kk*l*3),dim=c(l,kk,3))
  theta[1,,1]<-thetafirst
  theta[2,,1]<-oneiter(theta[1,,1])
  for(step in 2:(l-1)){
    theta[step+1,,1]=oneiter(theta[step,,1])
    theta[step-1,,2]=eps_1(theta[step-1,,1],theta[step,,1],theta[step+1,,1])
  }
  theta[1,,3]=eps_1(theta[1,,2],theta[2,,2],theta[3,,2])
  theta[2,,3]=eps_1(theta[2,,2],theta[3,,2],theta[4,,2])
  diff=metrics(theta[1,,3],theta[2,,3])
  while(diff>threshold){
    l=l+1
    theta=abind(theta,array(rep(0,kk*3),dim=c(1,kk,3)),along=1)
    theta[l,,1]=oneiter(theta[l-1,,1])
    theta[l-2,,2]=eps_1(theta[l-2,,1],theta[l-1,,1],theta[l,,1])
    theta[l-4,,3]=eps_1(theta[l-4,,2],theta[l-3,,2],theta[l-2,,2])
    diff=metrics(theta[l-4,,3],theta[l-5,,3])
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
      theta[step,,j]=theta[step+1,,j-2]+INV(theta[step+1,,j-1]-theta[step,,j-1])
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
        result_1=theta[l,,2]
      }else{
        theta[L,,k]=theta[L+1,,k-2]+INV(theta[L+1,,k-1]-theta[L,,k-1])
        if(k%%2==0){
          if(l%%2==0||(l%%2==1&&k<l)){
            diff_1=metrics(theta[L,,k],theta[L-1,,k])
            if(diff_1<diff){diff=diff_1;result_1=theta[L,,k]}
          }
        }
      }
    }
  }
  iternum[i,4]=l
  time[i,4]<-(proc.time()-d)[1]
}

par(mfrow = c(1, 2))
hist(iternum[,1],freq=TRUE,xlab="迭代次数",ylab="频数",main="em",col="red")
hist(iternum[,2],freq=TRUE,xlab="迭代次数",ylab="频数",main="eps",col="green")
hist(iternum[,3],freq=TRUE,xlab="迭代次数",ylab="频数",main="迭代delta2",col="blue")
hist(iternum[,4],freq=TRUE,xlab="迭代次数",ylab="频数",main="一般性epsilon",col="black")
hist(time[,1],freq=TRUE,xlab="迭代次数",ylab="频数",main="em",col="red")
hist(time[,2],freq=TRUE,xlab="迭代次数",ylab="频数",main="eps",col="green")
hist(time[,3],freq=TRUE,xlab="迭代次数",ylab="频数",main="迭代delta2",col="blue")
hist(time[,4],freq=TRUE,xlab="迭代次数",ylab="频数",main="一般性epsilon",col="black")


summary(iternum)
summary(time)






