

# Support for:
 # "NeuralSpatial: Statistics for Distinct neural mechanisms for heading retrieval and
 # context recognition in the hippocampus during spatial reorientation"
# author: "Manuel Miguel Ramos Álvarez and the Muzzio Lab"

  source("https://github.com/ManuMi68/StatMmRa/raw/main/RMmRaGen24.txt")
  source("https://github.com/ManuMi68/MuLaNa2/raw/main/Rallfun-v40vMM.txt")
  
  
  # Function chkPkg() is located in RMmRaGen.R
  chkPkg(c("data.table", "dplyr", "tidyr","tibble", "parallel", "foreach",
           "ggplot2","add2ggplot","introdataviz", "cowplot", "ggridges", "ggcorrplot",
           "ez","afex", "emmeans", "effectsize","scmamp",
           "rtf","mltools","apaTables",
           "nparLD", "nparcomp", "rankFD",
           "XNomial", "rstatix","EMT", "twosamples"))
  
  FxDir=getwd()
  DirPath<-"NeuroSpatial"
  if (!dir.exists(paste0(FxDir,"/",DirPath))) {dir.create(paste0(FxDir,"/",DirPath))}
  
  # Load all data
  load(gzcon(url("https://github.com/ManuMi68/MuLaNa2/raw/main/NeuroSpatialData/NeuroSpatial.RData")))
  
  # Color schemes for the 4 tasks
  # Tarea 1 (código 2 "Elevated Cliff")
  TskCol1<-c("green3","green4")
  TskCol2<-c("red","orange")
  TskCol3<-c("blue","dodgerblue")
  TskCol4<-c("gray28","gray50")
  
  TskCol1b<-c("red","orange","brown")
  
  
  # More specific functions, perhaps I will include them in RMMRA
  #define function to catch integer(0)
  integer0_test <- function(data) {
    if(identical(data, integer(0))) {
      return(0)
    }
    else {
      return(data)
    }}
  
  PosHocAut <-function(a0P, MxP, emP, OrdV,customP) {
    ExtMod<-integer0_test(grep("|",emP,fixed = T))
    em.simple.Global <- eval(parse(text = paste0("emmeans(a0P, ~ ", emP,")")))
    dtresAllGlobal<- setkeyv(as.data.table(em.simple.Global),OrdV)
    dtresAllGlobal[,`:=`(emmean=frmMM(emmean),SE=frmMM(SE,4),lower.CL=frmMM(lower.CL),upper.CL=frmMM(upper.CL))]
    ResPostH.Global<-pairs(em.simple.Global,adjust = "holm")
    ResPostH.GlobalOrto<-contrast(em.simple.Global, customP,adjust = "holm")
    ResPostH.Global.MM <-PosHoc.MM(em.simple.Global)     # Estimate Rom from emmeans
    em.Is<-em.r<-NA
    em.Is<- eval(parse(text = paste0("emmeans(a0P, ~ ", emP,",pbkrtest.limit = 18240)")))
    em.r<-  eval(parse(text = paste0("emmeans(MxP, ~ ", emP,",pbkrtest.limit = 18240)")))
    if (ExtMod==0) ResPosHAPA.Global <- ResAPA1.23(em.Is, em.r, MxP)
    if (!ExtMod==0) ResPosHAPA.Global <- ResAPA2.23(em.Is, em.r, MxP)
    #ResPosHAPA.Global <- ResAPA2.23(em.Is, em.r, MxP)
    Result<-list(
      Global=      em.simple.Global,
      Means=       dtresAllGlobal,
      PosHoc=      ResPostH.Global,
      PosHocPlan = ResPostH.GlobalOrto,
      PosHocMM =   ResPostH.Global.MM,
      PosHocAPA =  ResPosHAPA.Global
    )
  }
  
  ResAPA1.23<-function(em.simpleIsP, em.simple.rP, Mod.1rp) {
    pp1 =    pairs(em.simpleIsP,adjust = "holm")
    pp1Sum = summary(pp1)
    pp1.MM = PosHoc.MM(em.simpleIsP)
    pp1.r =  pairs(em.simple.rP,adjust = "holm")
    pp1b=    eff_size(em.simple.rP, sigma = sigma(Mod.1rp), edf = df.residual(Mod.1rp))
    pp1bSum= summary(pp1b)
    initEs=grep("estimate",names(data.table(pp1Sum)))-1
    
    # He sustituido pHolm por pRom
    # pp1.MM[1,pRom]
    # pp1Sum[1,"p.value"]
    ResAPA2<-list();#CtrEl=1
    for (CtrEl in 1:nrow(pp1Sum)) {
      ResAPA2[CtrEl]<- paste0(pp1Sum[CtrEl,1],": ",
                              "t(",formatC(round(pp1Sum[CtrEl,"df"],2),2,format="f"),") = ",
                              formatC(round(pp1Sum[CtrEl,"t.ratio"],2),2,format="f"), ", p = ",
                              # formatC(round(pp1Sum[CtrEl,"p.value"],4),4,format="f"), ", dPAIR = ",
                              pp1.MM[CtrEl,pRom], ", dPAIR = ",  
                              formatC(round(pp1bSum[CtrEl,"effect.size"],2),2,format="f"),
                              " (", InterpdCohen(pp1bSum[CtrEl,"effect.size"])," effect)",
                              ", CI 95% = [",
                              formatC(round(pp1bSum[CtrEl,"lower.CL"],2),2,format="f"), ", ",
                              formatC(round(pp1bSum[CtrEl,"upper.CL"],2),2,format="f"), "]"
      )}
    ResAPA2
  }
  
  ResAPA2.23<-function(em.simpleIsP, em.simple.rP, Mod.1rp) {
    pp1 =    pairs(em.simpleIsP,adjust = "holm")
    pp1Sum = summary(pp1)
    pp1.MM = PosHoc.MM(em.simpleIsP)
    pp1.r =  pairs(em.simple.rP,adjust = "holm")
    pp1b=    eff_size(em.simple.rP, sigma = sigma(Mod.1rp), edf = df.residual(Mod.1rp))
    pp1bSum= summary(pp1b)
    initEs=grep("estimate",names(data.table(pp1Sum)))-1
    
    ResAPA2<-list();#CtrEl=1
    for (CtrEl in 1:nrow(pp1Sum)) {
      ResAPA2[CtrEl]<- paste0(paste(unname(unlist((pp1Sum)[CtrEl,2:initEs])),collapse = " "),": ",pp1Sum[CtrEl,1],": ",
                              "t(",formatC(round(pp1Sum[CtrEl,"df"],2),2,format="f"),") = ",
                              formatC(round(pp1Sum[CtrEl,"t.ratio"],2),2,format="f"), ", p = ",
                              # formatC(round(pp1Sum[CtrEl,"p.value"],4),4,format="f"), ", dPAIR = ",
                              pp1.MM[CtrEl,pRom], ", dPAIR = ",  
                              formatC(round(pp1bSum[CtrEl,"effect.size"],2),2,format="f"),
                              " (", InterpdCohen(pp1bSum[CtrEl,"effect.size"])," effect)",
                              ", CI 95% = [",
                              formatC(round(pp1bSum[CtrEl,"lower.CL"],2),2,format="f"), ", ",
                              formatC(round(pp1bSum[CtrEl,"upper.CL"],2),2,format="f"), "]"
      )}
    ResAPA2
  }
  
  # Functions that computes Bayes Normandin-Ramos, adapted to the contrast of 0.25
  BayNormRam<-function(z,N,li=.25,ls=.9,Nulo=.25) {
    ct=1/(ls-li)
    Nu=ct*integrate(MiBer, z = z, N = N, lower = li,upper = ls)[[1]]
    De=MiBer(Nulo,z = z, N = N)
    Res=Nu/De
    Res}
  MiBer<-function(Theta,z,N) {(Theta^z) * ((1-Theta)^(N-z))}
  Productorio<-function(x) {Prodd=1;for (i in 1:length(x)) Prodd=Prodd*x[i];Prodd}
  
  BF1Grp <- function(DTp, colp=c("green3","green4"),logEs10=F) {
    ColPer<-c("cyan3","blue","red", "magenta2");
    if (logEs10) {LgEscala=log10(3);xlbl="Log10 (BF)"}
    if (!logEs10) {LgEscala=log(3); xlbl="Log (BF)"}
    P <- ecdf(DTp$logBF)
    Py = P(LgEscala)
    
    cdfbas<-ggplot (data=DTp, aes(x=logBF))
    cdfGrp<- cdfbas + stat_ecdf(linewidth=1,pad = T) +
      scale_colour_manual(values=colp) +
      theme_classic2() +
      theme(legend.title=element_blank(),
            legend.text=element_text(size=10),
            legend.position=c(0.2,0.8), 
            text = element_text(size=20),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            plot.caption = element_text(hjust = 0.5,size = 12)) +
      labs(y="Cumulative Proportion", x=xlbl,
           caption ="<----Favours chance/learning-----|-----Favours learning/chance---->") +
      xlim(-4,4)  +
      geom_vline(xintercept = LgEscala, linetype="dashed", linewidth=.25) +
      geom_vline(xintercept = -LgEscala, linetype="dashed", linewidth=.25) +
      geom_vline(xintercept = 0, linetype="dashed", linewidth=.25) +
      geom_hline(yintercept = .5, linetype="dashed", linewidth=.25) +
      annotate(geom="text", x=LgEscala, y=Py, label=percent(Py), color="red", hjust=0, vjust=0) +
      geom_point(aes(x=LgEscala,y=Py),color="red")
    cdfGrp
  }
  
  BFAllGrp1 <- function(DTpp, colp=c("green3","green4"),logEs10=F) {
    ColPer<-c("cyan3","blue","red", "magenta2");
    
    if (logEs10) {LgEscala=log10(3);xlbl="Log10 (BF)"}
    if (!logEs10) {LgEscala=log(3); xlbl="Log (BF)"}
    library(grid)
    
    JnGrp<-lapply(levels(DTpp$Day), function(xy) {
      DTppS<-DTpp[Day==xy]
      Py<-DTppS[, lapply(.SD, function(z) stats::ecdf(z)(LgEscala)), .SDcols = "logBF", by = .(Group) ]
      cdfbas<-ggplot (data=DTppS, aes(x=logBF, group =Group, color = Group, linetype = Group))
      cdfGrp<- cdfbas + stat_ecdf(linewidth=1,pad = T) +
        scale_colour_manual(values=colp) +
        theme_classic2() +
        theme(legend.title=element_blank(),
              legend.text=element_text(size=10),
              legend.position=c(0.2,0.8), 
              text = element_text(size=20),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.caption = element_text(hjust = 0.5,size = 12)) +
        labs(y="Cumulative Proportion", x=xlbl) +
        # labs(y="Cumulative Proportion", x=xlbl,
        #      caption ="<----Favours chance/learning-----|-----Favours learning/chance---->",cex.caption = 0.5) +
        xlim(-4,4)  +
        geom_vline(xintercept = LgEscala, linetype="dashed", linewidth=.25) +
        geom_vline(xintercept = -LgEscala, linetype="dashed", linewidth=.25) +
        geom_vline(xintercept = 0, linetype="dashed", linewidth=.25) +
        geom_hline(yintercept = .5, linetype="dashed", linewidth=.25) +
        annotate(geom="text", x=LgEscala, y=Py$logBF, label=percent(Py$logBF), color="gray50", hjust=0, vjust=0) +
        geom_point(data = data.frame(x=rep(LgEscala,length(levels(DTpp$Group))), Py),aes(x=x,y=logBF),color="gray50",size=2) +
        annotate("rect", xmin=-Inf, xmax=0, ymin=-Inf, ymax=1, alpha=0.2, fill="red") +
        annotate("rect", xmin=0, xmax=+Inf, ymin=-Inf, ymax=1, alpha=0.2, fill="green") +
        annotate(geom="text",x=4,y=.5,label="Favors Learning", cex=5, hjust=1, vjust=0) +
        annotate(geom="text",x=-4,y=.5,label="Favors Chance", cex=5, hjust=0, vjust=0) +
        coord_cartesian(ylim=c(0,1), clip = "off")
      cdfGrp
      
    })
    pltAll<-plot_grid(plotlist=JnGrp,nrow = 1, labels=levels(DTpp$Day), hjust=-2.5)
    names(JnGrp) <-levels(DTpp$Day)
    Res<-list(Session=JnGrp,All=pltAll )
    Res
  }
  
  BFAllGrp2 <- function(DTpp, colp=c("green3","green4"),logEs10=F) {
    ColPer<-c("cyan3","blue","red", "magenta2");
    if (logEs10) {LgEscala=log10(3);xlbl="Log10 (BF)"}
    if (!logEs10) {LgEscala=log(3); xlbl="Log (BF)"}
    
    JnGrp<-lapply(levels(DTpp$Group), function(xy) {
      DTppS<-DTpp[Group==xy]
      Py<-DTppS[, lapply(.SD, function(z) stats::ecdf(z)(LgEscala)), .SDcols = "logBF", by = .(Day) ]
      cdfbas<-ggplot (data=DTppS, aes(x=logBF, group =Day, color = Day, linetype = Day))
      cdfGrp<- cdfbas + stat_ecdf(linewidth=1,pad = T) +
        scale_colour_manual(values=colp) +
        theme_classic2() +
        theme(legend.title=element_blank(),
              legend.text=element_text(size=10),
              legend.position=c(0.2,0.8), 
              text = element_text(size=20),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.caption = element_text(hjust = 0.5,size = 12)) +
        labs(y="Cumulative Proportion", x=xlbl) +
        # labs(y="Cumulative Proportion", x=xlbl,
        #      caption ="<----Favours chance/learning-----|-----Favours learning/chance---->",cex.caption = 0.5) +
        xlim(-4,4)  +
        geom_vline(xintercept = LgEscala, linetype="dashed", linewidth=.25) +
        geom_vline(xintercept = -LgEscala, linetype="dashed", linewidth=.25) +
        geom_vline(xintercept = 0, linetype="dashed", linewidth=.25) +
        geom_hline(yintercept = .5, linetype="dashed", linewidth=.25) +
        annotate(geom="text", x=LgEscala, y=Py$logBF, label=percent(Py$logBF), color="gray50", hjust=0, vjust=0) +
        geom_point(data = data.frame(x=rep(LgEscala,length(levels(DTpp$Day))), Py),aes(x=x,y=logBF),color="gray50",size=2) +
        annotate("rect", xmin=-Inf, xmax=0, ymin=-Inf, ymax=1, alpha=0.2, fill="red") +
        annotate("rect", xmin=0, xmax=+Inf, ymin=-Inf, ymax=1, alpha=0.2, fill="green") +
        annotate(geom="text",x=4,y=.5,label="Favors Learning", cex=5, hjust=1, vjust=0) +
        annotate(geom="text",x=-4,y=.5,label="Favors Chance", cex=5, hjust=0, vjust=0) +
        coord_cartesian(ylim=c(0,1), clip = "off")
      cdfGrp 
      
    })
    pltAll<-plot_grid(plotlist=JnGrp,nrow = 1, labels=levels(DTpp$Group), hjust=-2.5)
    names(JnGrp) <-levels(DTpp$Group)
    Res<-list(Group=JnGrp,All=pltAll )
    Res
  }
  
  BFAllGrp3 <- function(DTpp, colp=c("green3","green4"),logEs10=F) {
    ColPer<-c("cyan3","blue","red", "magenta2");
    if (logEs10) {LgEscala=log10(3);xlbl="Log10 (BF)"}
    if (!logEs10) {LgEscala=log(3); xlbl="Log (BF)"}
    
    
    DTppS<-DTpp
    Py<-DTppS[, lapply(.SD, function(z) stats::ecdf(z)(LgEscala)), .SDcols = "logBF", by = .(Day) ]
    cdfbas<-ggplot (data=DTppS, aes(x=logBF, group =Day, color = Day, linetype = Day))
    cdfGrp<- cdfbas + stat_ecdf(linewidth=1,pad = T) +
      scale_colour_manual(values=colp) +
      theme_classic2() +
      theme(legend.title=element_blank(),
            legend.text=element_text(size=10),
            legend.position=c(0.2,0.8), 
            text = element_text(size=20),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            plot.caption = element_text(hjust = 0.5,size = 12)) +
      labs(y="Cumulative Proportion", x=xlbl) +
      # labs(y="Cumulative Proportion", x=xlbl,
      #      caption ="<----Favours chance/learning-----|-----Favours learning/chance---->",cex.caption = 0.5) +
      xlim(-4,4)  +
      geom_vline(xintercept = LgEscala, linetype="dashed", linewidth=.25) +
      geom_vline(xintercept = -LgEscala, linetype="dashed", linewidth=.25) +
      geom_vline(xintercept = 0, linetype="dashed", linewidth=.25) +
      geom_hline(yintercept = .5, linetype="dashed", linewidth=.25) +
      annotate(geom="text", x=LgEscala, y=Py$logBF, label=percent(Py$logBF), color="gray50", hjust=0, vjust=0) +
      geom_point(data = data.frame(x=rep(LgEscala,length(levels(DTpp$Day))), Py),aes(x=x,y=logBF),color="gray50",size=2) +
      annotate("rect", xmin=-Inf, xmax=0, ymin=-Inf, ymax=1, alpha=0.2, fill="red") +
      annotate("rect", xmin=0, xmax=+Inf, ymin=-Inf, ymax=1, alpha=0.2, fill="green") +
      annotate(geom="text",x=4,y=.5,label="Favors Learning", cex=5, hjust=1, vjust=0) +
      annotate(geom="text",x=-4,y=.5,label="Favors Chance", cex=5, hjust=0, vjust=0) +
      coord_cartesian(ylim=c(0,1), clip = "off")
    cdfGrp 
    
    
    
  }
  
  AcumBFLog.23<-function(DTp,Gs,Ds,logEs10=T,wG=T,wS=T) {
    ValsBF2=seq(-2, 2, by=1.0)
    
    DTp<-DTp[Group==Gs&Day==Ds]
    Global  <- Productorio(DTp$BF)
    ResBFg1 <- DTp$logBF
    P <- ecdf(DTp$logBF)
    
    if (logEs10) {          
      ValsBF=sort(c(seq(-2,2,by=.5),-log10(3),log10(3),-.01));
      GlobalLog<-log10(Global)
    }
    if (!logEs10) {
      ValsBF=sort(c(seq(-2,2,by=.5),-log(3),log(3),-.01));
      GlobalLog<-log(Global)
    }        
    
    if (wG&wS) laEtiq=paste0(Gs," ",Ds)
    if (wG&!wS) laEtiq=Gs
    if (!wG&wS) laEtiq=Ds
    
    BFAc <-list(
      BF =       DTp$BF,
      InterBF =  InterpBF(DTp$BF),
      LogBF =    DTp$logBF,	
      DT =       data.table(Group=rep(laEtiq,length(ResBFg1)),BF=ResBFg1),
      Pecdf =    P,
      empcdf =   empirical_cdf(ResBFg1, ubounds=ValsBF),
      empcdf2 =  empirical_cdf(ResBFg1, ubounds=ValsBF2),
      NegPerc =  percent(P(0)),
      BFGlobal = paste0("BF All= ",round(Global,3), "; Log BF All= ",round(GlobalLog,3), "; ",InterpBF(Global), " Evidence"),
      BFAll =    round(Global,3),
      LogBFAll = round(GlobalLog,3),
      InterBFAll = InterpBF(Global),
      Grp <- BF1Grp(DTp = DTp,logEs10 = logEs10)
    )
    BFAc
  }
  
  AcumBFLog.23.NoCx<-function(DTp,Ds,logEs10=T,wG=T,wS=T) {
    ValsBF2=seq(-2, 2, by=1.0)
    
    DTp<-DTp[Day==Ds]
    Global  <- Productorio(DTp$BF)
    ResBFg1 <- DTp$logBF
    P <- ecdf(DTp$logBF)
    
    if (logEs10) {          
      ValsBF=sort(c(seq(-2,2,by=.5),-log10(3),log10(3),-.01));
      GlobalLog<-log10(Global)
    }
    if (!logEs10) {
      ValsBF=sort(c(seq(-2,2,by=.5),-log(3),log(3),-.01));
      GlobalLog<-log(Global)
    }        
    
    laEtiq=Ds
    
    BFAc <-list(
      BF =       DTp$BF,
      InterBF =  InterpBF(DTp$BF),
      LogBF =    DTp$logBF,	
      DT =       data.table(BF=ResBFg1),
      Pecdf =    P,
      empcdf =   empirical_cdf(ResBFg1, ubounds=ValsBF),
      empcdf2 =  empirical_cdf(ResBFg1, ubounds=ValsBF2),
      NegPerc =  percent(P(0)),
      BFGlobal = paste0("BF All= ",round(Global,3), "; Log BF All= ",round(GlobalLog,3), "; ",InterpBF(Global), " Evidence"),
      BFAll =    round(Global,3),
      LogBFAll = round(GlobalLog,3),
      InterBFAll = InterpBF(Global),
      Grp <- BF1Grp(DTp = DTp,logEs10 = logEs10)
    )
    BFAc
  }
  
  com.Glob <- function(DTp) {
    res<- foreach(i=levels(DTp$Group),.combine='c') %:%
      foreach(j=levels(DTp$Day)) %dopar% { 
        AcumBFLog.23(DTp, i, j, F, T, T) }
    names(res) <-(foreach(i=levels(DTp$Group),.combine='c') %:%
                    foreach(j=levels(DTp$Day), .combine='c') %dopar% {
                      paste0(i, ".", j)
                    })
    res
  }
  
  com.Glob.NoCx <- function(DTp) {
    foreach(j=levels(DTp$Day)) %dopar% { 
      AcumBFLog.23.NoCx(DTp, j, F, T, T) }
  }
  
  com.Glob.APA <- function(DTp) {
    foreach(i=levels(DTp$Group),.combine='c') %:%
      foreach(j=levels(DTp$Day), .combine='c') %dopar% {
        paste0(i," ",j,": ",AcumBFLog.23(DTp, i, j, F, T, T)$BFGlobal)}
  }
  
  com.Glob.APA.NoCx <- function(DTp) {
    foreach(j=levels(DTp$Day), .combine='c') %dopar% {
      paste0(j,": ",AcumBFLog.23.NoCx(DTp, j, F, T, T)$BFGlobal)}
  }
  
  com.ecdf<-function(DTp) {
    Resecdf<-list()
    setkey(DTp,Group,Day)
    lv1=length(levels(DTp$Group));lv2=length(levels(DTp$Day))
    Resecdf[[1]]<-  unlist(lapply(1:lv1, function (x) lapply(1:lv2, function (y) 
      DTp[.(levels(DTp$Group)[x],levels(DTp$Day)[y],with=F),stats::ecdf(logBF)])),recursive = F)
    
    names(Resecdf[[1]]) <-(foreach(i=levels(DTp$Group),.combine='c') %:%
                             foreach(j=levels(DTp$Day), .combine='c') %dopar% {
                               paste0(i, ".", j)
                             })
    
    # Resecdf[[2]]<-  unlist(lapply(1:lv1, function (x) lapply(1:lv2, function (y) 
    #  DTp[.(levels(DTp$Group)[x],levels(DTp$Day)[y],with=F),empirical_cdf(logBF, ubounds=ValsBF)])),
    #  recursive = F)
    # Resecdf[[3]]<-  unlist(lapply(1:lv1, function (x) lapply(1:lv2, function (y) 
    #  DTp[.(levels(DTp$Group)[x],levels(DTp$Day)[y],with=F),empirical_cdf(logBF, ubounds=ValsBF2)])),
    #  recursive = F)
    Resecdf[[2]] <-foreach(i=1:lv1,.combine='cbind') %:%
      foreach(j=1:lv2, .combine='cbind') %dopar% {
        sapply (DTp[.(levels(DTp$Group)[i],levels(DTp$Day)[j],with=F),
                    empirical_cdf(logBF, ubounds=ValsBF)]$CDF, percent)
      } %>% as.data.frame %>%
      `colnames<-`(foreach(i=levels(DTp$Group),.combine='c') %:%
                     foreach(j=levels(DTp$Day), .combine='c') %dopar% {
                       paste0(i, ".", j)
                     }) %>%
      mutate(UpperBound= round(ValsBF,2), .before=colnames(.)[1])
    Resecdf[[3]] <-foreach(i=1:lv1,.combine='cbind') %:%
      foreach(j=1:lv2, .combine='cbind') %dopar% {
        sapply (DTp[.(levels(DTp$Group)[i],levels(DTp$Day)[j],with=F),
                    empirical_cdf(logBF, ubounds=ValsBF2)]$CDF, percent)
      } %>% as.data.frame %>%
      `colnames<-`(foreach(i=levels(DTp$Group),.combine='c') %:%
                     foreach(j=levels(DTp$Day), .combine='c') %dopar% {
                       paste0(i, ".", j)
                     }) %>%
      mutate(UpperBound= round(ValsBF2,2), .before=colnames(.)[1])
    Resecdf
  }
  
  com.ecdf.NoCx<-function(DTp) {
    Resecdf<-list()
    setkey(DTp,Day)
    lv1=1;lv2=length(levels(DTp$Day))
    Resecdf[[1]]<-  unlist(lapply(1:lv1, function (x) lapply(1:lv2, function (y) 
      DTp[.(levels(DTp$Day)[y],with=F),stats::ecdf(logBF)])),recursive = F)
    names(Resecdf[[1]]) <-levels(DTp$Day)
    Resecdf[[2]] <-foreach(i=1:lv1,.combine='cbind') %:%
      foreach(j=1:lv2, .combine='cbind') %dopar% {
        sapply (DTp[.(levels(DTp$Day)[j],with=F),
                    empirical_cdf(logBF, ubounds=ValsBF)]$CDF, percent)
      } %>% as.data.frame %>%
      `colnames<-`(
        foreach(j=levels(DTp$Day), .combine='c') %dopar% {
          paste0(j)
        }) %>%
      mutate(UpperBound= round(ValsBF,2), .before=colnames(.)[1])
    Resecdf[[3]] <-foreach(i=1:lv1,.combine='cbind') %:%
      foreach(j=1:lv2, .combine='cbind') %dopar% {
        sapply (DTp[.(levels(DTp$Day)[j],with=F),
                    empirical_cdf(logBF, ubounds=ValsBF2)]$CDF, percent)
      } %>% as.data.frame %>%
      `colnames<-`(
        foreach(j=levels(DTp$Day), .combine='c') %dopar% {
          paste0(j)
        }) %>%
      mutate(UpperBound= round(ValsBF2,2), .before=colnames(.)[1])
    Resecdf
  }
  
  MakeBF<-function(DTppp, Nmpp) {
    Global.i<- DTppp %>% 
      .[,.(BF=Productorio(BF)),by=list(Group,Day)] %>%
      .[,InterBF:=InterpBF(BF)] %>%
      .[,logBF:=log(BF)]  %>%
      .[,BF:=frmMM(BF,4)] %>%
      .[,logBF:=frmMM(logBF,4)] %>%
      #dplyr::rename(., BF=V1) %>%
      data.table()
    ecdf.i<-com.ecdf(DTppp)
    Detall.i <- com.Glob(DTppp)
    GlobalAPA.i <- com.Glob.APA(DTppp)
    GrpGroup.i<-BFAllGrp1(DTppp,TskCol2)
    GrpDay.i<-  BFAllGrp2(DTppp, TskCol1b)
    GrpJn.i<-plot_grid(plotlist=lapply (Detall.i, "[[", 13), 
                       labels = c(t(outer(levels(DTppp$Group), levels(DTppp$Day), FUN=paste,sep=" "))))
    
    write.csv2(Global.i,paste0(Nmpp," Global.csv"))
    zz<-file(paste0(Nmpp," Global.txt"),"w")
    sink(zz);
    print(GlobalAPA.i);
    cat("\n");
    print(ecdf.i[[2]]);
    cat("\n"); cat("\n");
    print(ecdf.i[[3]]);
    cat("\n");cat("\n");cat("\n");
    print(Detall.i)
    sink(); close(zz)
    for (i in 1:length(levels(DTppp$Group))) {ggsave(paste0(Nmpp," CumDis_",levels(DTppp$Group)[i], ".pdf"),GrpDay.i[[1]][[i]])}
    for (i in 1:length(levels(DTppp$Day)))   {ggsave(paste0(Nmpp," CumDis_",levels(DTppp$Day)[i], ".pdf"),GrpGroup.i[[1]][[i]])}
    ggsave(paste0(Nmpp," CumDis_CxtBySess.pdf"),GrpGroup.i[[2]],width = 12,height = 6)
    ggsave(paste0(Nmpp," CumDis_SessByCxt.pdf"),GrpDay.i[[2]],width = 8,height = 6)
    ggsave(paste0(Nmpp," CumDis_All.pdf"),GrpJn.i,width = 20,height = 14)
    
    Res <-list(
      Global= Global.i,
      ecdfMM = ecdf.i,
      Detall = Detall.i,
      GlobalAPA = GlobalAPA.i,
      GrpGroup = GrpGroup.i,
      GrpDay = GrpDay.i,
      GrpJn = GrpJn.i
    )
    Res
  }
  
  MakeBFNoCx<-function(DTppp, Nmpp) {
    Global.i<- DTppp %>% 
      .[,.(BF=Productorio(BF)),by=list(Day)] %>%
      .[,InterBF:=InterpBF(BF)] %>%
      .[,logBF:=log(BF)]  %>%
      .[,BF:=frmMM(BF,4)] %>%
      .[,logBF:=frmMM(logBF,4)] %>%
      #dplyr::rename(., BF=V1) %>%
      data.table()
    ecdf.i<-com.ecdf.NoCx(DTppp)
    Detall.i <- com.Glob.NoCx(DTppp)
    names(Detall.i)<-levels(DTppp$Day)
    GlobalAPA.i <- com.Glob.APA.NoCx(DTppp)
    # GrpGroup.i<-BFAllGrp1(DTppp,TskCol2)
    GrpDay.i<-  BFAllGrp3(DTppp, TskCol1b)
    GrpJn.i<-plot_grid(plotlist=lapply (Detall.i, "[[", 13), 
                       labels = levels(DTppp$Day))
    
    write.csv2(Global.i,paste0(Nmpp," Global.csv"))
    zz<-file(paste0(Nmpp," Global.txt"),"w")
    sink(zz);
    print(GlobalAPA.i);
    cat("\n");
    print(ecdf.i[[2]]);
    cat("\n"); cat("\n");
    print(ecdf.i[[3]]);
    cat("\n");cat("\n");cat("\n");
    print(Detall.i)
    sink(); close(zz)
    ggsave(paste0(Nmpp," CumDis_SessByCxt.pdf"),GrpDay.i)
    ggsave(paste0(Nmpp," CumDis_All.pdf"),GrpJn.i,width = 20,height = 14)
    
    Res <-list(
      Global= Global.i,
      ecdfMM = ecdf.i,
      Detall = Detall.i,
      GlobalAPA = GlobalAPA.i,
      # GrpGroup = GrpGroup.i,
      GrpDay = GrpDay.i,
      GrpJn = GrpJn.i
    )
    Res
  }
  
  stats_MM <- function(x) {
    nadt <-0
    summ<-summary(x)
    if (length(summ)>6) nadt= summ[[7]]
    n = length(x)-nadt
    min <- frmMM(min(x, na.rm = TRUE),2)
    max <- frmMM(max(x, na.rm = TRUE),2)
    mean <- frmMM(mean(x, na.rm = TRUE),2)
    med <- frmMM(median(x, na.rm = TRUE))
    trimm = frmMM(mean(x, na.rm = TRUE,trim=.2),2)
    Q1 <- frmMM(summ[[2]],2)
    Q3 <- frmMM(summ[[5]],2)
    Iqr <- frmMM(summ[[5]]-summ[[2]],2)
    sd <- frmMM(sd(x, na.rm =TRUE),2)
    se <- frmMM(sd(x, na.rm =TRUE)/sqrt(n),2)
    mad <- frmMM(mad(x, na.rm = TRUE),2)
    summary <- data.table(n=n, 'NAs' = nadt, Min = min, '1st Qu'=Q1, 
                          Median=med, Mean = mean, 'Trimmed(20%)' = trimm,
                          '3rd Qu'=Q3, Max = max, 
                          SD = sd, SEM = se, IQR = Iqr, MAD = mad)
    summary
  }
  
  ResCtrl1p <- function(DTp, DTRes, a0, PosH.4w) {
    # Results
    cat("Data Structure:\n");
    str(DTp) # Data File
    cat("---------------------\n");
    print(kableTabl(DTp,"Data", ""))
    cat("---------------------\n");
    print(kableTabl(ezPrecis(DTp),"Design Structure", ""))
    cat("---------------------\n");
    print(kableTabl(DTRes,"Descriptive", ""))
    cat("Omnibus AOV:\n")
    print(a0)
    cat("---------------------\n");
    cat("Pos Hoc Simple Effects:\n")
    print(PosH.4w$PosHocMM)
    #kableTabl(PosH.4w$PosHocMM,"Pos Hoc Simple Effects", "")
    cat("---------------------\n");
    print(kableTabl(PosH.4w$Means,"Descrptive & CI-95%", ""))
    cat("---------------------\n");
    cat("Graphics have been preloaded to avoid page overload\n (see the corresponding section for details).:\n");
    # plot(Grp) # Exploratory Analysis
  }
  
  InterP2 <- function(x) {
    dplyr::case_when(
      x <= .001   ~ "***",
      x <= .01    ~ "** ",
      x <= .05    ~ "*  ",
      x <  .07     ~ "#  ",
      .default =    "   ")
  }
  
  MkMtrSig <-function(lv1=2,lv2=4) {
    NmLev=lapply (1:lv1, function(ii) paste0("C",ii,".",1:lv2))
    MtrL=matrix("",nrow = lv2+1,ncol=lv2+1)
    rownames(MtrL) <-c(NmLev[[1]],"Marg")
    colnames(MtrL) <-c(NmLev[[2]],"Marg")
    MtrL
  }
  
  LimpMtr<-function(Mtrx) {Mtrx[is.na(Mtrx)]<-""; Mtrx[Mtrx=="ns"]<-""; Mtrx}
  
  LimpMtr2 <- function(Mtrx) {Mtrx<- gsub("\\s", "", Mtrx); Mtrx[is.na(Mtrx)]<-""; Mtrx[Mtrx=="ns"]<-""; Mtrx[Mtrx=="NA"]<-""; Mtrx}
  
  CreateDT.Means3 <- function(Filep , ip='1',jp="across",kp='FS') {
    DTCom=NA
    DTCom <- Filep %>% filter(Day==jp,Context==kp,CellType==ip) %>% droplevels() %>% 
      select(-c(Day:CellType,Context )) %>%
      data.table() 
    
    DTCom2 <- DTCom %>%
      .[, data.table(table(Subject,bmr)),] %>%
      .[,Ntt:=sum(N),   by=(Subject)] %>% 
      .[,Perc:=N/sum(N),by=(Subject)] %>%
      .[,BMROrder:=paste0("C",rank(-Perc, ties.method= "first")),by=.(Subject)] %>%
      data.table() %>%
      setkey(., Subject) %>%
      filter(BMROrder %in% c("C1","C2")) %>% droplevels() %>% 
      mutate_at(c("BMROrder","bmr"), factor) %>%
      data.table() %>%
      .[,Inter:=interaction(BMROrder,bmr,lex.order = T)] %>% setkey(.,Inter) %>% 
      mutate(Inter=factor(Inter,)) %>% data.table()%>%
      .[,Condition:=paste(Inter,collapse = "_"), by=.(Subject)] %>%  
      mutate(Condition=factor(Condition,levels=LosLevLim)) %>%
      data.table() 
    
    DTAnimalFus=data.table(left_join(DTCom2,DTCom[,.(Subject,animal)], by="Subject",multiple="first"))
    
    # Overall estimates
    DTComAll <- DTAnimalFus %>%
      .[, data.table(table(Condition,animal)),] %>%
      .[,Ntt:=sum(N), by=.(animal)] %>%
      .[,Perc:=N/Ntt] %>%
      data.table()
    
    
    pWilcox <-c(do.call("rbind",lapply(1:12, function(ll) wilcox.test(DTComAll[Condition==LosLevLim[[ll]]]$Perc,
                                                                      mu = 1/12, alternative = "greater")$p.value)))
    Individual.a=data.table(Condition=LosLevLim , pWilcox)
    
    
    Fus2=full_join(DTComAll[,.(Mean=mean(Perc)), by =.(Condition)],Individual.a,by= "Condition")
    
    Resf<-as.data.table(Fus2) %>%
      .[,pw.Sig:= InterP2(pWilcox)] %>% 
      data.table()
    
    # Marginal estimations
    DTCom3aa = DTAnimalFus %>% filter(BMROrder %in% c("C1")) %>% droplevels() %>% 
      mutate_at(c("BMROrder","bmr"), factor) %>% data.table() %>%
      .[,data.table(table(bmr, animal))] %>%
      mutate(bmr=factor(bmr)) %>% as.data.table() %>%
      .[,Ntt:=sum(N), by=.(animal)] %>%
      .[,Perc:=N/Ntt] %>%
      data.table()
    
    DTCom3bb = DTAnimalFus %>% filter(BMROrder %in% c("C2")) %>% droplevels() %>% 
      mutate_at(c("BMROrder","bmr"), factor) %>% data.table() %>%
      .[,data.table(table(bmr, animal))] %>%
      mutate(bmr=factor(bmr)) %>% as.data.table() %>%
      .[,Ntt:=sum(N), by=.(animal)] %>%
      .[,Perc:=N/Ntt] %>%
      data.table()
    
    
    Individual.aaa <- Individual.bbb <-NA
    
    pWilcox.aaa <- c(do.call("rbind",lapply(1:4, function(ll) wilcox.test(DTCom3aa[bmr==levels(bmr)[[ll]]]$Perc,
                                                                          mu = 1/4, alternative = "greater")$p.value)))
    Individual.aaa=data.table(bmr=levels(DTCom3aa$bmr), pWilcox=pWilcox.aaa)
    Individual.aaa=full_join(DTCom3aa[,.(Mean=mean(Perc)), by =.(bmr)],Individual.aaa,by= "bmr")
    
    pWilcox.bbb <- c(do.call("rbind",lapply(1:4, function(ll) wilcox.test(DTCom3bb[bmr==levels(bmr)[[ll]]]$Perc,
                                                                          mu = 1/4, alternative = "greater")$p.value)))
    Individual.bbb=data.table(bmr=levels(DTCom3bb$bmr), pWilcox=pWilcox.bbb)
    Individual.bbb=full_join(DTCom3bb[,.(Mean=mean(Perc)), by =.(bmr)],Individual.bbb,by= "bmr")
    
    Individual.aaa=cbind(Condition="C1",Individual.aaa)
    Individual.bbb=cbind(Condition="C2",Individual.bbb)
    
    Fus111=rbind(Individual.aaa, Individual.bbb)
    
    Resf2<- as.data.table(Fus111) %>%
      .[,pw.Sig:= InterP2(pWilcox)] %>% 
      data.table()
    
    Resfinal <-list(
      MedGlobal <- Resf,
      MedMargC1 <-as.data.table(Resf2[Condition=="C1"]),
      MedMargC2 <-as.data.table(Resf2[Condition=="C2"])
    )
    
    Resfinal
  }
