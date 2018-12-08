library(dplyr)
library(ggplot2)
library(lme4)
library(stargazer)
library(Hmisc)
library(retimes)
library(dotwhisker)
require(RColorBrewer)
library("ggpubr")

## LOAD AND CLEAN DATA
Book1=read.csv("/Users/jasondsc/Desktop/McGill/EEG_alertness/Behav_EEG_alert/good/alertness_cleaned_behav_2018.txt",sep = "\t" )
Book1=Book1[,1:19]
Book1=Book1[Book1$Tgt.Type!=0,]
Book1=Book1[Book1$Disc.Key.Press==70 |Book1$Disc.Key.Press==74 , ]
Book1=Book1[Book1$Subj.Key.Press==70 |Book1$Subj.Key.Press==74 , ]
Book1$Subject.Code[Book1$Subject.Code=="051_AG_2"]="051_AG"
#Book1=Book1[!(Book1$Subject.Code %in% c("014_AR", "041_PN", "043_JB", "048_JS", "051_AG", "052_SM", "054_WX_2", "055_FH", "057_EK")),]
Book1=Book1[Book1$Disc.Key.RT<2.5,]
Book1=Book1[Book1$Disc.Key.RT>0.150,]
Book1$TM.SOA[Book1$TM.SOA==2]=0
Book1$TM.SOA[Book1$TM.SOA==4]=1
Book2= Book1[Book1$Disc.Key.Accurac==1,]
book3=Book1

## CALCULATE D_PRIME AND C for GLM
df = Book1 %>% group_by(Subject.Code, TM.SOA, CuePresent,Tgt.Type ) %>%
  summarise(Avg=mean(Disc.Key.Accurac), N=n())

df$Avg[df$Tgt.Type==2] = 1-df$Avg[df$Tgt.Type==2]
df$Avg[df$Avg == 0] = 1/(2*df$N[df$Avg == 0])
df$Avg[df$Avg == 1] = ((2*df$N[df$Avg == 1])-1)/(2*df$N[df$Avg == 1])

df$Z.Avg= qnorm(df$Avg, 0, 1)

temp = df$Z.Avg[df$Tgt.Type==1] -  df$Z.Avg[df$Tgt.Type==2]
criterion = -0.5*(df$Z.Avg[df$Tgt.Type==1] + df$Z.Avg[df$Tgt.Type==2])

df2 = df[which(df$Tgt.Type==1),]
df2$d.prime = temp
df2$C= criterion
df2$Disc.Key.Expecte = c()
df2$Z.Avg = c()

df2$TM.SOA[df2$TM.SOA==2]=0
df2$TM.SOA[df2$TM.SOA==4]=1

df2$CuePresent=as.factor(df2$CuePresent)
df2$TM.SOA=as.factor(df2$TM.SOA)

## SET VARIABLES OF INTREST TO FACTORS

Book1$CuePresent=as.factor(Book1$CuePresent)
Book1$Subj.Report=as.factor(Book1$Subj.Report)
Book1$TM.SOA=as.factor(Book1$TM.SOA)
Book2$CuePresent=as.factor(Book2$CuePresent)
Book2$TM.SOA=as.factor(Book2$TM.SOA)

## CALCULATE EXGAUSSIAN PARAMETERS
Exgauss= Book2 %>% group_by(Subject.Code, TM.SOA, CuePresent) %>%
  summarise(MU=mexgauss(Disc.Key.RT)[[1]], SIGMA=mexgauss(Disc.Key.RT)[[2]], TAU=mexgauss(Disc.Key.RT)[[3]])

Exgauss= Book2 %>% group_by(Subject.Code, TM.SOA, CuePresent) %>%
  summarise(MU=timefit(Disc.Key.RT)@par[[1]], SIGMA=timefit(Disc.Key.RT)@par[[2]], TAU=timefit(Disc.Key.RT)@par[[3]])

Exgauss$CuePresent=as.factor(Exgauss$CuePresent)
Exgauss$TM.SOA=as.factor(Exgauss$TM.SOA)

## CALCULATE D PRIME AND EXGAUSS FOR CORRELATIONS
cue_Exgauss= Book2 %>% group_by(Subject.Code, CuePresent) %>%
  summarise(MU=timefit(Disc.Key.RT)@par[[1]], SIGMA=timefit(Disc.Key.RT)@par[[2]], TAU=timefit(Disc.Key.RT)@par[[3]])

cue_Exgauss= Book2 %>% group_by(Subject.Code, CuePresent, Subj.Report) %>%
  summarise(MU=timefit(Disc.Key.RT)@par[[1]], SIGMA=timefit(Disc.Key.RT)@par[[2]], TAU=timefit(Disc.Key.RT)@par[[3]])

df3 = Book1 %>% group_by(Subject.Code, CuePresent,Tgt.Type ) %>%
  summarise(Avg=mean(Disc.Key.Accurac), N=n())

df3$Avg[df3$Tgt.Type==2] = 1-df3$Avg[df3$Tgt.Type==2]
df3$Avg[df3$Avg == 0] = 1/(2*df3$N[df3$Avg == 0])
df3$Avg[df3$Avg == 1] = ((2*df3$N[df3$Avg == 1])-1)/(2*df3$N[df3$Avg == 1])

df3$Z.Avg= qnorm(df3$Avg, 0, 1)

cue_dprime = df3$Z.Avg[df3$Tgt.Type==1] -  df3$Z.Avg[df3$Tgt.Type==2]
cue_criterion = -0.5*(df3$Z.Avg[df3$Tgt.Type==1] + df3$Z.Avg[df3$Tgt.Type==2])

df4 = df3[which(df3$Tgt.Type==1),]
df4$d.prime = cue_dprime
df4$C= cue_criterion

df5 = Book1 %>% group_by(Subject.Code, CuePresent, Subj.Report,Tgt.Type ) %>%
  summarise(Avg=mean(Disc.Key.Accurac), N=n())

df5$Avg[df5$Tgt.Type==2] = 1-df5$Avg[df5$Tgt.Type==2]
df5$Avg[df5$Avg == 0] = 1/(2*df5$N[df5$Avg == 0])
df5$Avg[df5$Avg == 1] = ((2*df5$N[df5$Avg == 1])-1)/(2*df5$N[df5$Avg == 1])

df5$Z.Avg= qnorm(df5$Avg, 0, 1)

cue_dprime = df5$Z.Avg[df5$Tgt.Type==1] -  df5$Z.Avg[df5$Tgt.Type==2]
cue_criterion = -0.5*(df5$Z.Avg[df5$Tgt.Type==1] + df5$Z.Avg[df5$Tgt.Type==2])

df6 = df5[which(df5$Tgt.Type==1),]
df6$d.prime = cue_dprime
df6$C= cue_criterion


################### ERPs BRAIN SCORES CORRELATIONS ###########################################

cueaccsub = Book1 %>% group_by(Subject.Code, CuePresent, Subj.Report ) %>%
  summarise(Avg=mean(Disc.Key.Accurac), N=n())

erp=read.csv("~/Desktop/pls_analysis_eeg_alertness/ERP_spectral_brainscroes.csv", sep = ",", header = FALSE)
cueaccsub = Book1 %>% group_by(Subject.Code, CuePresent, Subj.Report ) %>%
  summarise(Avg=mean(Disc.Key.Accurac), N=n())

plot(cueaccsub$Avg[(cueaccsub$CuePresent==1 & cueaccsub$Subj.Report==0)], erp$V1[1:30])
plot(cueaccsub$Avg[(cueaccsub$CuePresent==0 & cueaccsub$Subj.Report==0)], erp$V1[31:60])
plot(cueaccsub$Avg[(cueaccsub$CuePresent==1 & cueaccsub$Subj.Report==1)], erp$V1[61:90])
plot(cueaccsub$Avg[(cueaccsub$CuePresent==0 & cueaccsub$Subj.Report==1)], erp$V1[91:120])

cor.test(cueaccsub$Avg[(cueaccsub$CuePresent==1 & cueaccsub$Subj.Report==0)], erp$V1[1:30])
cor.test(cueaccsub$Avg[(cueaccsub$CuePresent==0 & cueaccsub$Subj.Report==0)], erp$V1[31:60])
cor.test(cueaccsub$Avg[(cueaccsub$CuePresent==1 & cueaccsub$Subj.Report==1)], erp$V1[61:90])
cor.test(cueaccsub$Avg[(cueaccsub$CuePresent==0 & cueaccsub$Subj.Report==1)], erp$V1[91:120])

cor.test(df6$d.prime[(df6$CuePresent==1 & df6$Subj.Report==0)], erp$V1[1:30])
cor.test(df6$d.prime[(df6$CuePresent==0 & df6$Subj.Report==0)], erp$V1[31:60])
cor.test(df6$d.prime[(df6$CuePresent==1 & df6$Subj.Report==1)], erp$V1[61:90])
cor.test(df6$d.prime[(df6$CuePresent==0 & df6$Subj.Report==1)], erp$V1[91:120])

data1=data.frame(df6$d.prime[(df6$CuePresent==0 & df6$Subj.Report==1)], erp$V1[91:120])
colnames(data1)= c("C", "Brain")

data2=data.frame(df6$d.prime[(df6$CuePresent==1 & df6$Subj.Report==1)], erp$V1[91:120])
colnames(data2)= c("C", "Brain")

ggscatter(data1, x= "Brain", y="C",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Brain Scores", ylab = "D prime")

ggscatter(data2, x= "Brain", y="C",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Brain Scores", ylab = "D prime")


data1=data.frame(df6$d.prime[(df6$CuePresent==0 & df6$Subj.Report==1)], erp$V1[91:120])
colnames(data1)= c("C", "Brain")

data2=data.frame(df6$d.prime[(df6$CuePresent==1 & df6$Subj.Report==1)], erp$V1[61:90])
colnames(data2)= c("C", "Brain")

ggscatter(data1, x= "Brain", y="C",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Brain Scores", ylab = "D prime")

ggscatter(data2, x= "Brain", y="C",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Brain Scores", ylab = "D prime")


################### SPECTRAL BRAIN SCORES CORRELATIONS ###########################################

spect=read.csv("~/Desktop/pls_analysis_eeg_alertness/cueing_spectral_brainscroes.csv", sep = ",", header = FALSE)
spect$log_1=log(spect$V1*10^10)

cueSuj = book3 %>% group_by(Subject.Code, CuePresent ) %>%
  summarise(Avg=mean(Subj.Report), N=n())

cueacc = book3 %>% group_by(Subject.Code, CuePresent ) %>%
  summarise(Avg=mean(Disc.Key.Accurac), N=n())

cueaccsub = Book1 %>% group_by(Subject.Code, CuePresent, Subj.Report ) %>%
  summarise(Avg=mean(Disc.Key.Accurac), N=n())

cue_Exgauss= Book2 %>% group_by(Subject.Code, CuePresent) %>%
  summarise(MU=timefit(Disc.Key.RT)@par[[1]], SIGMA=timefit(Disc.Key.RT)@par[[2]], TAU=timefit(Disc.Key.RT)@par[[3]])
cor.test(spect[1:30,1][-28], df4$d.prime[df4$CuePresent==1][-28])
cor.test(spect[31:60,1][-28], df4$d.prime[df4$CuePresent==0][-28])
cor.test(spect[1:30,1][-28], cueSuj$Avg[cueSuj$CuePresent==1][-28])
cor.test(spect[31:60,1][-28], cueSuj$Avg[cueSuj$CuePresent==0][-28])
cor.test(spect[1:30,1][-28], cue_Exgauss$MU[cue_Exgauss$CuePresent==1][-28])
cor.test(spect[31:60,1][-28], cue_Exgauss$MU[cue_Exgauss$CuePresent==0][-28])


################### GLMER ###########################################

lm0= lmer(d.prime ~ (1 | Subject.Code), df2)
lm1=lmer(d.prime ~ (1 | Subject.Code) + CuePresent, df2)
lm2=lmer(d.prime ~ (1 | Subject.Code) + TM.SOA, df2)
lm3=lmer(d.prime ~ (1 | Subject.Code) + CuePresent + TM.SOA, df2)
lm4=lmer(d.prime ~ (1 | Subject.Code) + CuePresent+TM.SOA+CuePresent:TM.SOA, df2)

anova(lm0, lm1, lm3,lm4)
summary(lm3)
confint(lm3)

lm0= lmer(C ~ (1 | Subject.Code), df2)
lm1=lmer(C ~ (1 | Subject.Code) + CuePresent, df2)
lm3=lmer(C ~ (1 | Subject.Code) + CuePresent + TM.SOA, df2)
lm4=lmer(C ~ (1 | Subject.Code) + CuePresent+TM.SOA+CuePresent:TM.SOA, df2)

anova(lm0, lm1, lm3,lm4)

lm0= lmer(Disc.Key.RT ~ (1 | Subject.Code), Book2)
lm1=lmer(Disc.Key.RT ~ (1 | Subject.Code) + CuePresent, Book2)
lm3=lmer(Disc.Key.RT ~ (1 | Subject.Code) + CuePresent + TM.SOA, Book2)
lm4=lmer(Disc.Key.RT ~ (1 | Subject.Code) + CuePresent+TM.SOA+CuePresent:TM.SOA, Book2)

anova(lm0, lm1,lm3,lm4)

lm0= lmer(MU ~ (1 | Subject.Code), Exgauss)
lm1=lmer(MU ~ (1 | Subject.Code) + CuePresent, Exgauss)
lm3=lmer(MU ~ (1 | Subject.Code) + CuePresent + TM.SOA, Exgauss)
lm4=lmer(MU ~ (1 | Subject.Code) + CuePresent+TM.SOA+CuePresent:TM.SOA, Exgauss)

anova(lm0, lm1,lm3,lm4)
summary(lm3)
confint(lm3)

lm0= lmer(SIGMA ~ (1 | Subject.Code), Exgauss)
lm1=lmer(SIGMA ~ (1 | Subject.Code) + CuePresent, Exgauss)
lm3=lmer(SIGMA ~ (1 | Subject.Code) + CuePresent + TM.SOA, Exgauss)
lm4=lmer(SIGMA ~ (1 | Subject.Code) + CuePresent+TM.SOA+CuePresent:TM.SOA, Exgauss)

anova(lm0, lm1,lm3,lm4)

lm0= lmer(TAU ~ (1 | Subject.Code), Exgauss)
lm1=lmer(TAU ~ (1 | Subject.Code) + CuePresent, Exgauss)
lm3=lmer(TAU ~ (1 | Subject.Code) + CuePresent + TM.SOA, Exgauss)
lm4=lmer(TAU ~ (1 | Subject.Code) + CuePresent+TM.SOA+CuePresent:TM.SOA, Exgauss)

anova(lm0, lm1,lm3,lm4)
summary(lm3)
confint(lm3)

lm0= glmer(Subj.Report ~ (1 | Subject.Code), Book1, family = "binomial")
lm1=glmer(Subj.Report ~ (1 | Subject.Code) + CuePresent, Book1, family = "binomial")
lm3=glmer(Subj.Report ~ (1 | Subject.Code) + CuePresent + TM.SOA, Book1, family = "binomial")
lm4=glmer(Subj.Report ~ (1 | Subject.Code) + CuePresent+TM.SOA+CuePresent:TM.SOA, Book1, family = "binomial")

anova(lm0, lm1,lm3,lm4)
summary(lm3)
confint(lm3)

lm0= glmer(Disc.Key.Accurac ~ (1 | Subject.Code), Book1, family = "binomial")
lm1=glmer(Disc.Key.Accurac ~ (1 | Subject.Code) + CuePresent, Book1, family = "binomial")
lm3=glmer(Disc.Key.Accurac ~ (1 | Subject.Code) + CuePresent + TM.SOA, Book1, family = "binomial")
lm4=glmer(Disc.Key.Accurac ~ (1 | Subject.Code) + CuePresent+TM.SOA+ Subj.Report, Book1, family = "binomial")
lm5=glmer(Disc.Key.Accurac ~ (1 | Subject.Code) + CuePresent+TM.SOA+CuePresent:TM.SOA +Subj.Report, Book1, family = "binomial")
lm6=glmer(Disc.Key.Accurac ~ (1 | Subject.Code) + CuePresent+TM.SOA+ Subj.Report +Subj.Report:CuePresent +CuePresent:TM.SOA , Book1, family = "binomial")
lm7=glmer(Disc.Key.Accurac ~ (1 | Subject.Code) + CuePresent+TM.SOA+ Subj.Report +Subj.Report:CuePresent+Subj.Report:TM.SOA +CuePresent:TM.SOA , Book1, family = "binomial")
lm8=glmer(Disc.Key.Accurac ~ (1 | Subject.Code) + CuePresent+TM.SOA+ Subj.Report +Subj.Report:CuePresent+Subj.Report:TM.SOA+Subj.Report:TM.SOA:CuePresent +CuePresent:TM.SOA , Book1, family = "binomial")
anova(lm0, lm1,lm3,lm4,lm5,lm6,lm7,lm8)
summary(lm4)
ll=confint(lm4)

############################ GGPLOT FIGURES OF GLM RESULTS ###########################################

mean11=df2 %>% group_by(CuePresent, TM.SOA) %>% summarise(avg=mean(d.prime), N=n(), sterr=sd(d.prime)/sqrt(N), CI_low=mean_cl_boot(d.prime)[[2]],CI_high=mean_cl_boot(d.prime)[[3]])

ggplot(mean11, aes(CuePresent, avg, group=TM.SOA, color=TM.SOA, size)) + geom_point(size=5)+ geom_line(size=2) +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), mean11, width=0.2, size=2) + theme_minimal() +labs(y= "d prime", x ="") +  
  scale_color_brewer(palette="Pastel1", name  ="Masking \n Latency", breaks=c("0", "1"), labels=c("26 ms", "52 ms")) +
  scale_x_discrete(breaks=c("0","1"), labels=c("Not Cued", "Cued")) + theme(legend.justification=c(1,1), legend.position=c(1,1))

mean11=Book1 %>% group_by(CuePresent, TM.SOA, Subj.Report) %>% summarise(avg=mean(Disc.Key.Accurac), N=n(), sterr=sd(Disc.Key.Accurac)/sqrt(N), CI_low=mean_cl_boot(Disc.Key.Accurac)[[2]],CI_high=mean_cl_boot(Disc.Key.Accurac)[[3]])
mean11$TM.SOA = plyr::revalue(mean11$TM.SOA, c("0"="26 ms", "1"="52 ms"))
mean11$Subj.Report = plyr::revalue(mean11$Subj.Report, c("0"="Not Seen", "1"="Seen"))
mean11 %>% ggplot(aes(CuePresent, avg, fill = Subj.Report)) + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0.1, size=0.75, position = position_dodge(0.9))+ 
  scale_fill_brewer(palette="Pastel1") + theme_minimal() + labs(y= "Discrimination \n Accuracy", x ="", fill= "Subjective \n Report")  + scale_x_discrete(breaks=c("0","1"), labels=c("Not Cued", " Cued")) + facet_wrap(~TM.SOA ) +
  coord_cartesian(ylim=c(0.5, 0.9))

mean11=Exgauss %>% group_by(CuePresent, TM.SOA) %>% summarise(avg=mean(MU), N=n(), sterr=sd(MU)/sqrt(N), CI_low=mean_cl_boot(MU)[[2]],CI_high=mean_cl_boot(MU)[[3]])
mean11$TM.SOA=as.factor(mean11$TM.SOA)
mean11$CuePresent=as.factor(mean11$CuePresent)
mean11$TM.SOA = plyr::revalue(mean11$TM.SOA, c("0"="26 ms", "1"="52 ms"))
ggplot(mean11, aes(CuePresent, avg, group=TM.SOA, fill=TM.SOA, size)) + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0.1, size=0.75, position = position_dodge(0.9))+ 
  scale_fill_brewer(palette="Pastel1") + theme_minimal() + labs(y= "Mu", x ="", fill= "Masking \n Latency")  + scale_x_discrete(breaks=c("0","1"), labels=c("Not Cued", " Cued")) +
  coord_cartesian(ylim=c(0.3, 0.6)) + theme(legend.justification=c(1,1), legend.position="none")


mean11=Exgauss %>% group_by(CuePresent, TM.SOA) %>% summarise(avg=mean(TAU), N=n(), sterr=sd(TAU)/sqrt(N), CI_low=mean_cl_boot(TAU)[[2]],CI_high=mean_cl_boot(TAU)[[3]])
mean11$TM.SOA=as.factor(mean11$TM.SOA)
mean11$CuePresent=as.factor(mean11$CuePresent)
mean11$TM.SOA = plyr::revalue(mean11$TM.SOA, c("0"="26 ms", "1"="52 ms"))
ggplot(mean11, aes(CuePresent, avg, group=TM.SOA, fill=TM.SOA, size)) + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0.1, size=0.75, position = position_dodge(0.9))+ 
  scale_fill_brewer(palette="Pastel1") + theme_minimal() + labs(y= "Tau", x ="", fill= "Masking \n Latency")  + scale_x_discrete(breaks=c("0","1"), labels=c("Not Cued", " Cued")) +
  coord_cartesian(ylim=c(0, 0.4)) +theme(legend.justification=c(1,1), legend.position="none")


mean11=df2 %>% group_by(CuePresent, TM.SOA) %>% summarise(avg=mean(d.prime), N=n(), sterr=sd(d.prime)/sqrt(N), CI_low=mean_cl_boot(d.prime)[[2]],CI_high=mean_cl_boot(d.prime)[[3]])
mean11$TM.SOA=as.factor(mean11$TM.SOA)
mean11$CuePresent=as.factor(mean11$CuePresent)
mean11$TM.SOA = plyr::revalue(mean11$TM.SOA, c("0"="26 ms", "1"="52 ms"))
ggplot(mean11, aes(CuePresent, avg, group=TM.SOA, fill=TM.SOA, size)) + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0.1, size=0.75, position = position_dodge(0.9))+ 
  scale_fill_brewer(palette="Pastel1") + theme_minimal() + labs(y= "d prime", x ="", fill= "d prime")  + scale_x_discrete(breaks=c("0","1"), labels=c("Not Cued", " Cued")) +
  coord_cartesian(ylim=c(0, 3)) +theme(legend.justification=c(1,1), legend.position="none")

mean11=book3 %>% group_by(CuePresent, TM.SOA) %>% summarise(avg=mean(Subj.Report), N=n(), sterr=sd(Subj.Report)/sqrt(N), CI_low=mean_cl_boot(Subj.Report)[[2]],CI_high=mean_cl_boot(Subj.Report)[[3]])
mean11$CuePresent=as.factor(mean11$CuePresent)
mean11$TM.SOA=as.factor(mean11$TM.SOA)
mean11$TM.SOA = plyr::revalue(mean11$TM.SOA, c("0"="26 ms", "1"="52 ms"))
ggplot(mean11, aes(CuePresent, avg, group=TM.SOA, fill=TM.SOA, size)) + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0.1, size=0.75, position = position_dodge(0.9))+ 
  scale_fill_brewer(palette="Pastel1") + theme_minimal() + labs(y= "Subjective \n Report", x ="", fill= "d prime")  + scale_x_discrete(breaks=c("0","1"), labels=c("Not Cued", " Cued")) +
  coord_cartesian(ylim=c(0, 1)) +theme(legend.justification=c(1,1), legend.position="none")


########################################### MAKE CONTRAST PLOTS ########################################

con1= read.csv("~/Desktop/contrast2.csv", header = FALSE)
con1$group=as.factor(c("Cued","Not Cued","Cued","Not Cued"))
colnames(con1)[2:3]=c("CI_low", "CI_high")

ggplot(con1, aes(1:4, V1, group=group, fill=group, size)) + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=CI_low, ymax=CI_high),width=0.1, size=0.75, position = position_dodge(0.9) )+
  scale_fill_brewer(palette="Pastel1") + theme_minimal() + labs(y= "Optimal Contrast", x ="", fill= "")  +
  coord_cartesian(ylim=c(-0.45*10^-4, 0.45*10^-4)) + scale_x_discrete(labels=c("1.5" = "Not Seen"))

ggplot(con1, aes(1:4, V2, group=group, fill=group, size)) + geom_bar(stat="identity", position="dodge") +
  scale_fill_brewer(palette="Pastel1") + theme_minimal() + labs(y= "Optimal Contrast", x ="", fill= "")  +
  coord_cartesian(ylim=c(-0.6, 0.6)) + scale_x_discrete(labels=c("1.5" = "Not Seen"))

con1= read.csv("~/Desktop/pls_analysis_eeg_alertness/behav_spect_contrast.csv", header = FALSE)
con1$group=as.factor(c("Cued","Not Cued"))

ggplot(con1, aes(1:2, V1, group=group, fill=group, size)) + geom_bar(stat="identity", position="dodge") +
  scale_fill_brewer(palette="Pastel1") + theme_minimal() + labs(y= "Optimal Contrast", x ="", fill= "")  +
  coord_cartesian(ylim=c(-0.65, 0.65)) + scale_x_discrete(labels=c("1.5" = "Not Seen"))

