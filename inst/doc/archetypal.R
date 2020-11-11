## ----setup, include=FALSE-----------------------------------------------------
library(archetypal)
library(rmarkdown)
knitr::opts_chunk$set(echo = TRUE)
options(max.width = 1000)
options(max.print = 100000)

## ----2D, echo=TRUE------------------------------------------------------------
library(archetypal)
p1=c(1,2);p2=c(3,5);p3=c(7,3) 
dp=rbind(p1,p2,p3);dp
set.seed(9102)
pts=t(sapply(1:100, function(i,dp){
  cc=runif(3)
  cc=cc/sum(cc)
  colSums(dp*cc)
},dp))
df=data.frame(pts)
colnames(df)=c("x","y")
head(df)

## ----plot2D, echo=FALSE,out.width='75%', fig.align='center',fig.width=10, fig.height=6----

plot(dp,pch=3,cex=2,xlab='x',ylab='y',xlim=c(min(dp[,1])-0.25,max(dp[,1])+0.25),ylim=c(min(dp[,2])-0.25,max(dp[,2])+0.25))
points(df,col='blue',pch=19,cex=0.7)
polygon(rbind(dp,dp[1,]),lty=2)

## ----load2D, echo=TRUE--------------------------------------------------------
# data("wd2")
# df=wd2

## ----run2D, echo=TRUE---------------------------------------------------------
aa = archetypal(df = df, kappas = 3, verbose = TRUE, rseed = 9102, save_history = TRUE)

## ----out2D, echo=TRUE---------------------------------------------------------
names(aa)


## ----archplot, echo=FALSE,out.width='75%', fig.align='center',fig.width=10, fig.height=6----
plot(dp,pch=3,cex=2,xlab='x',ylab='y',xlim=c(min(dp[,1])-0.25,max(dp[,1])+0.25),ylim=c(min(dp[,2])-0.25,max(dp[,2])+0.25))
points(df,col='blue',pch=19,cex=0.7)
polygon(rbind(dp,dp[1,]),lty=2)
archs=data.frame(aa$BY)
points(archs,col='blue',pch=15,cex=2)
polygon(rbind(archs,archs[1,]),col=rgb(0, 0, 1,0.5)) 

## ----sse_conv, echo=FALSE,out.width='75%', fig.align='center',fig.width=10, fig.height=6----
vsse=aa$run_results$SSE
plot(vsse,xlab="Iteration",ylab="SSE",pch=19, col="blue",type="b")
grid()

## ----checkB, echo=TRUE--------------------------------------------------------
BB=aa$B
yy=check_Bmatrix(B = BB, chvertices = NULL, verbose = TRUE)
# yy$used_rows
# yy$used_weights

## ----ch2d, echo=TRUE----------------------------------------------------------
ch=chull(df)
ch
df[ch,]

## ----checkBCH, echo=TRUE------------------------------------------------------
yy$used_rows
unlist(yy$used_rows)%in%ch

## ----Barchplot, echo=FALSE,out.width='75%', fig.align='center',fig.width=10, fig.height=6----
plot(dp,pch=3,cex=2,xlab='x',ylab='y',xlim=c(min(dp[,1])-0.25,max(dp[,1])+0.25),ylim=c(min(dp[,2])-0.25,max(dp[,2])+0.25))
points(df,col='blue',pch=19,cex=0.7)
polygon(rbind(dp,dp[1,]),lty=2)
archs=data.frame(aa$BY)
points(archs,col='blue',pch=15,cex=2)
pp=lapply(yy$used_rows,function(x,df){points(df[x,],col='red',type='b',pch=19);lines(df[x,],col='red',lwd=2)},df)

## ----run2Dinitial, echo=TRUE--------------------------------------------------
aa2=archetypal(df=df,kappas = 3,initialrows =  c(34,62,86), verbose = TRUE,rseed=9102,save_history = TRUE)
yy2=check_Bmatrix(aa2$B,verbose = TRUE)

## ----plot3D, echo=TRUE,out.width='100%', fig.align='center',fig.width=17, fig.height=9----
library(plot3D)
#
p1=c(3,0,0);p2=c(0,5,0);p3=c(3,5,7);p4=c(0,0,0);
dp=data.frame(rbind(p1,p2,p3,p4));dp=dp[chull(dp),];colnames(dp)=c("x","y","z")
set.seed(9102)
df=data.frame(t(sapply(1:100, function(i,dp){
  cc=runif(4)
  cc=cc/sum(cc)
  colSums(dp*cc)
},dp)))
colnames(df)=c("x","y","z")
scatter3D(x=dp$x,y=dp$y,z=dp$z,colvar=NULL,lwd = 2, d = 3,xlab='x',ylab='y',zlab='z',theta=120,phi=15,
          main = "Generators and Data Points", bty ="g",ticktype = "detailed",col='black',pch=10,cex=2.5)
points3D(x=df$x,y=df$y,z=df$z,col='blue',add=T,pch=19)

## ----load3D, echo=TRUE--------------------------------------------------------
# data("wd3")
# df=wd3

## ----run3D, echo=TRUE---------------------------------------------------------
aa3 = archetypal(df = df, kappas = 4, verbose = TRUE, rseed = 9102, save_history = TRUE)
yy3 = check_Bmatrix(aa3$B)

## ----run3DuseBetas, echo=TRUE-------------------------------------------------
irows=yy3$leading_rows
aa4 = archetypal(df = df, kappas = 4, initialrows = irows, verbose = TRUE, rseed = 9102, save_history = TRUE)
yy4 = check_Bmatrix(aa4$B)

## ----plot3DBetas, echo=FALSE,out.width='100%', fig.align='center',fig.width=17, fig.height=9----
scatter3D(x=dp$x,y=dp$y,z=dp$z,colvar=NULL,lwd = 2, d = 3,xlab='x',ylab='y',zlab='z',theta=120,phi=15,
          main = "Archetypes and Used Points", bty ="g",ticktype = "detailed",col='black',pch=10,cex=2.5)
points3D(x=df$x,y=df$y,z=df$z,col='blue',add=TRUE,pch=19)
archs3=data.frame(aa4$BY)
points3D(archs3$x,archs3$y,archs3$z,col='blue',add=TRUE,pch=15,cex=2.5)
pp3=lapply(yy4$used_rows,function(x,df){
  dh=df[x,]
  points3D(x=dh$x,y=dh$y,z=dh$z,col='red',add=TRUE,pch=19,cex=1.5)
  if(length(x)!=1){segments3D(x0=dh$x[1],y0=dh$y[1],z0=dh$z[1],x1=dh$x[2],y1=dh$y[2],z1=dh$z[2],col='red',add=TRUE,lwd=3) }
  },df)

## ---- ch3,echo=TRUE-----------------------------------------------------------
ch=unique(do.call(c,as.list(geometry::convhulln(df,'Fx'))))
ch

## ----find1, echo=TRUE---------------------------------------------------------
yy1 = find_outmost_projected_convexhull_points(df, kappas = 4)
yy1$outmost
yy1$outmostall
yy1$outmostall%in%ch

## ----find2, echo=TRUE---------------------------------------------------------
yy2 = find_outmost_convexhull_points(df, kappas = 4)
yy2$outmost
yy2$outmostall
yy2$outmostall%in%ch

## ----find3, echo=TRUE---------------------------------------------------------
# yy3 = find_outmost_partitioned_convexhull_points(df, kappas = 4, nworkers = 10)
# yy3$outmost
# yy3$outmostall
# yy3$outmostall%in%ch
# 1 . 2 . 3 . 4 . 5 . 6 . 7 . 8 . 9 . 10 .   
# Time difference of 2.769091 secs
# [1] 84  3
# [1] 61 64 82 67
# [1] 61 64 82 67
# [1] TRUE TRUE TRUE TRUE

## ----find4, echo=TRUE---------------------------------------------------------
# yy4 = find_furthestsum_points(df, kappas = 4, nfurthest = 100, nworkers = 10, sortrows = TRUE)
# yy4$outmost
# yy4$outmostall
# yy4$outmostall%in%ch
# [1] 56 61 64 67
# [1] 56 61 64 67
# [1] TRUE TRUE TRUE TRUE

## ----find5, echo=TRUE---------------------------------------------------------
yy5 = find_outmost_points(df, kappas = 4)
yy5$outmost
yy5$outmostall
yy5$outmostall%in%ch

