# author: "Manuel Miguel Ramos Álvarez and the Muzzio Lab" 2024

AdjMod.23.f<-function(LasMedP,TyLbl=NULL,xLb="X",yLb="Y",xmaxp=10,yminp=0,PIp=NULL,
                      LineRate=TRUE,pointXMed="p",pointXMedLb="OnAx",RootG=FALSE,MaxG=FALSE,InflG=FALSE) {
  chkPkg(c("shape","RootsExtremaInflections"))
  decilecolors=c("#7F7FCE", "#7F7FF7", "#88A9F9", "#93D2FB", "#A0FCFE",
                 "#BDFDD7", "#DDFEB3", "#FFFF91", "#F8D68B", "#F3AE86" )
  colorsdef=rep(decilecolors, each=6)
  #ResF<-list()
  NmMod="RegAs1"
  
  XMin=head(LasMedP$X,1)
  XMax=tail(LasMedP$X,1)
  YMax=LasMedP[LasMedP$X==XMin,]$Y;
  YMin=LasMedP[LasMedP$X==XMax,]$Y;
  RatePrv=-(YMax-YMin)/(XMax-XMin);
  # xPred <- seq(XMin, XMax, length = 1000);
  
  fit<-NA
  try(fit<-nlsAsym2(YMin,YMax,RatePrv,LasMedP))
  
  xPred <- seq(0, XMax, length = 1000);
  PredMod=predict(fit, data.frame(X = xPred))
  plot(xPred, PredMod,lty=1,lwd=3,xlab=xLb,ylab=yLb,type="l",xlim=c(0,xmaxp),ylim=c(yminp, 1),frame.plot=FALSE,cex.lab=1.25)
  if (!is.null(TyLbl)) {mtext(do.call(expression, TyLbl),side=3,line=(length(TyLbl)-1):0)}
  
  tetha1=a=Asymptote= Asym=AsymStim=coefficients(fit)[[1]]
  theta2=b=Origin=    R0=R0Stim=coefficients(fit)[[2]]
  theta3=LogRate=     lrc=lrcStim=coefficients(fit)[[3]]
  xMed=log(2)/(exp(lrcStim))
  PredxMed=predict(fit, data.frame(X = xMed))
  Rate=c=             exp(lrcStim)
  
  # Find root, plot results, print Taylor coefficients and rho estimation:
  bR <-rootxi(LasMedP$X,LasMedP$Y,1,length(LasMedP$X),5,5,plots=F); # bR$froot[2]
  # Find extreme, plot results, print Taylor coefficients and rho estimation:
  cR <-extremexi(LasMedP$X,LasMedP$Y,1,length(LasMedP$X),5,5,plots=F); # c$fextr[2]
  # Find inflection point, plot results, print Taylor coefficients and rho estimation:
  dR <-inflexi(LasMedP$X,LasMedP$Y,1,length(LasMedP$X),5,5,plots=F);   # d$finfl[2]
  
  if (RootG | MaxG | InflG) {
    PredModRoot<-rootxiGraph(LasMedP$X,LasMedP$Y,1,length(LasMedP$X),5,5,TRUE)
    lines(LasMedP$X, PredModRoot, lty=1,lwd=2,cex.lab=1.25,col="red")
    legend(2,yminp+.2, col = c("red"), lty = 1, 
           lwd = 1, legend = c(paste0("Taylor fit (Polynomial Order 5)")), bty = "n", cex = 0.7)
    legend(2,yminp+.25, col = c("black"), lty = 1, 
           lwd = 3, legend = bquote(.("Asymptotic Regression")*~(RSE==.(round(SEMod(fit),4)))), bty = "n", cex = 0.7)
  } else{
    with(LasMedP, lines(X,Y,type="b",pch=16));
  }
  
  # Las Medias al final
  with(LasMedP, points(X,Y,pch=16, cex = 1.45));
  with(LasMedP, points(X,Y,pch=16, col=decilecolors,cex = 1.25));
  
  #log(2)
  if (LineRate) {
    newx <- 1
    pred0 <- data.frame(x=newx, y=AsymRegMM(newx,AsymStim,R0Stim,c))
    pred1 <- data.frame(x=newx, y=DrvAsymReg(newx,AsymStim,R0Stim,c))
    yint <- pred0$y - (pred1$y*newx)
    xint <- -yint/pred1$y
    lines(xPred, yint + pred1$y*xPred, lty=2,lwd=1) # tangent (1st deriv. of spline at newx)
    #points(xint, 0, col=3, pch=19) # x intercept
    LTang=lm(yint + pred1$y*xPred ~ xPred)$coefficients
    ratio=5/1
    Angle = atan(LTang[2] * ratio) * (180 / pi) # Yo estimé 70
    text(2, .75,"Max Rate at x = 1",srt=Angle,adj=c(0.75,0),cex=.90) 
  }
  
  if (pointXMed=="p") points(xMed, PredxMed[[1]],pch=16,cex=1.5,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5))
  if (pointXMed=="X") points(xMed, PredxMed[[1]],pch=4, cex=1.5,col=rgb(red = 0, green = 1, blue = 1))
  if (pointXMed %in% c("p", "X") & pointXMedLb =="OnAx") {
    lines(c(0,xMed),c(PredxMed[[1]],PredxMed[[1]]),lty=2,lwd=1)
    lines(c(xMed,xMed),c(PredxMed[[1]],-1),lty=2,lwd=1)
    # text(xMed+3, PredxMed,expression(X[0.5]%~~%.(xMed)))
    # text(xMed+1.5, PredxMed,bquote(X[0.5]==.(ceiling(xMed))),cex=1.5)
    text(xMed+.5, PredxMed[[1]],
         quote({x[0.5] %~~%1}), 
         adj = c(0, 0.5),
         cex=1)
    abline(h=AsymStim,lty=2,lwd=1)
    abline(h=R0Stim,lty=2,lwd=1)
    axis(side=2, at = PredxMed[[1]], labels= "",las=1, cex=.75)
    axis(side=1, at = xMed, labels= "",las=1, cex=.75)
    par(xpd=TRUE)
    text(x = -1.3, y= .1,labels=expression(frac(theta[1]+theta[2],2)),adj=c(.5,0),cex=.75)
    text(x = xMed, y= -1.4,labels=expression(frac(log(2), e^{theta[3]})==frac(log(2), Rate)),cex=.75)
    Arrows(-1,.15,-.4,PredxMed[[1]], arr.length = .25, arr.type="triangle", arr.adj = 1)
    Arrows(1,-1.35,1,-1.1, arr.length = .25, arr.type="triangle", arr.adj = 1 )
    par(xpd=FALSE)
    axis(side=2, at = R0Stim, labels= expression(theta[2]), pos=-.5, las=1)
    axis(side=2, at = AsymStim, labels= expression(theta[1]), pos=-.5, las=1)
    # axis(2, at = PredxMed[[1]], labels= expression(frac(theta[1]+theta[2],2)), 
    #                pos=-.25, las=1,cex=2)
  }
  
  if (pointXMed %in% c("p", "X") & pointXMedLb =="OnGrp") {
    lines(c(0,xMed),c(PredxMed[[1]],PredxMed[[1]]),lty=2,lwd=1)
    lines(c(xMed,xMed),c(PredxMed[[1]],-1),lty=2,lwd=1)
    text(xMed+.5, PredxMed[[1]],
         quote({f(x[0.5]) == frac(theta[1]+theta[2],2)}~{} %=>% {}~~
                 {x[0.5] == frac(log(2), e^{theta[3]})}~ {} %=>% {}~~ {x[0.5] %~~%1}), 
         adj = c(0, 0.5),
         cex=1)
    text(x = 0.5, y= PredxMed[[1]],labels=expression(frac(theta[1]+theta[2],2)),cex=.75)
  }
  
  if(!is.null(PIp)) {
    abline(v=PIp,lty=3,lwd=1)
    text(PIp, 1,"IP",cex=1)
  }
  
  if (RootG) {
    abline(v=bR$froot[2],lty=3,lwd=1)
    text(bR$froot[2],1, paste0("Root = ",bR$froot[2]),cex=1)
  }
  if (MaxG) {
    abline(v=cR$fextr[2],lty=3,lwd=1)
    text(cR$fextr[2],1, paste0("Extreme = ",cR$fextr[2]),cex=1)
  }
  if (InflG){
    abline(v=dR$finfl[2],lty=3,lwd=1)
    text(dR$finfl[2],1, paste0("Inflection = ",dR$finfl[2]),cex=1)
  }
  
  NNN.1 <-  AdjParms.LM(fit)
  for(ii in (1:3)) NNN.1[[ii]] <- round(NNN.1[[ii]],2)
  ResM.1= data.table(GenTypeM="E.Mod",ClasType="AsymRegres",
                     Model="SSasymp Orig",as.data.table(NNN.1))
  ResM.2= paste0(paste0(names(unlist(NNN.1[1:3])),"= ",  unlist(NNN.1[1:3]),collapse="; "),"; ",
                 NNN.1[[4]])
  
  ResF <-list(
    Type=NmMod,
    XPred=xPred,
    Model=fit,
    Rsq=round(1 - var(residuals(fit))/var(LasMedP$Y),4),
    SE=SEMod(fit),
    Summ = summary(fit),
    Par= c(coefficients(fit),Rate=Rate,XHalf=xMed,y_Xhalf=(AsymStim+R0Stim)/2),
    forAPA1=ResM.1,
    forAPA2=ResM.2,
    BondAj=AdjParms(fit,LasMedP$Y),
    root = list(an=bR$an, froot =bR$froot),
    extr = list(an =cR$an, extr =cR$fextr),
    inflexi = list(an =dR$an, finfl =dR$finfl)
  )
  ResF
}
