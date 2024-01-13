packages <- c("survival","ggplot2","ggfortify")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
require("survival")
require("ggplot2")
require("ggfortify")


tempos=c(1, 2, 3, 3, 3, 5, 5, 16,16,16, 16, 16, 16, 16, 16,
1, 1, 1, 1, 4, 5, 7, 8, 10, 10, 12, 16, 16, 16)
censura=c(0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
1,1,1,0,0,1,1,1,1,0,0,0,0,0)
grupo=c("Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle",
"Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide")

ekm <- survfit(Surv(tempos,censura)~grupo)
summary(ekm)
autoplot(ekm)

s5=(1-3/14)*(1-1/9)
s5

vs5=((1-3/14)*(1-1/9))^2*(3/(14*(14-3))+1/(9*(9-1)))
ds5=sqrt(vs5)
ds5

ekm <- survfit(Surv(tempos,censura)~grupo,conf.type="plain")
summary(ekm)

ekm <- survfit(Surv(tempos,censura)~grupo,conf.type="log-log") #Kalbfleish e Prentice(1980)
summary(ekm)

ekm <- survfit(Surv(tempos,censura)~grupo,conf.type="log")
summary(ekm)

ekm <- survfit(coxph(Surv(tempos[grupo=="Esteróide"],censura[grupo=="Esteróide"])~1,method="breslow")) #Nelson-Aalen
summary(ekm)
racum=-log(ekm$surv)
racum

s10tv=(1-3/(14-1/2*2))*(1-3/(9-1/2*0)) #tabela de vida
s10tv

####----------Estimação de quantidades básicas---------####
#exemplo de hepatite
packages <- c("survival","ggplot2","ggfortify")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
require("survival")
require("ggplot2")
require("ggfortify")

tempos=c(1, 2, 3, 3, 3, 5, 5, 16,16,16, 16, 16, 16, 16, 16,
1, 1, 1, 1, 4, 5, 7, 8, 10, 10, 12, 16, 16, 16)
censura=c(0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
1,1,1,0,0,1,1,1,1,0,0,0,0,0)
grupo=c("Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle",
"Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide")

censE=censura[grupo=="Esteróide"]
tempE=tempos[grupo=="Esteróide"]
ekmE=survfit(Surv(tempE,censE)~1)
summary(ekmE)

summary(ekmE)$surv[summary(ekmE)$time==7]

#sobrevivência: interpolação
s6=(6-5)*(summary(ekmE)$surv[summary(ekmE)$time==7]-summary(ekmE)$surv[summary(ekmE)$time==5])/(7-5)+summary(ekmE)$surv[summary(ekmE)$time==5]
s6

#tempo mediano de vida: interpolação
t50=(10-8)*(0.5-summary(ekmE)$surv[summary(ekmE)$time==8])/(summary(ekmE)$surv[summary(ekmE)$time==10]-summary(ekmE)$surv[summary(ekmE)$time==8])+8
t50

#tempo médio de vida
t=tempE[censE==1]
tj=c(0,as.numeric(levels(as.factor(t))))
surv=c(1,as.numeric(levels(as.factor(ekmE$surv))))
surv=sort(surv,decreasing=T)
k=length(tj)-1
prod=matrix(0,k,1)
	for(j in 1:k){
		prod[j]=(tj[j+1]-tj[j])*surv[j]
		}
tm=sum(prod)
tm

####----------Reincidência de tumor sólido---------####
packages <- c("survival","ggplot2","ggfortify")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
require("survival")
require("ggplot2")
require("ggfortify")

tempos=c(3,4,5.7,6.5,6.5,8.4,10,10,12,15)
censura=c(1,0,0,1,1,0,1,0,1,1)

ekm <- survfit(Surv(tempos,censura)~1,conf.type="plain")
summary(ekm)
autoplot(ekm)

#tempo mediano de vida: interpolação
t50=(10-6.5)*(0.5-summary(ekm)$surv[summary(ekm)$time==6.5])/(summary(ekm)$surv[summary(ekm)$time==10]-summary(ekm)$surv[summary(ekm)$time==6.5])+6.5
t50

#tempo médio de vida
t=tempos[censura==1]
tj=c(0,as.numeric(levels(as.factor(t))))
surv=c(1,as.numeric(levels(as.factor(ekm$surv))))
surv=sort(surv,decreasing=T)
k=length(tj)-1
prod=matrix(0,k,1)
	for(j in 1:k){
		prod[j]=(tj[j+1]-tj[j])*surv[j]
		}
tm=sum(prod)
tm

#vida média residual
vmr10=((15-12)*0.241+(12-10)*0.482)/0.482
vmr10

#-------------------Comparação de curvas de sobrevivência----------------#
#LOGRANK
packages <- c("survival","ggplot2","ggfortify")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
require("survival")
require("ggplot2")
require("ggfortify")

tempos=c(1, 2, 3, 3, 3, 5, 5, 16,16,16, 16, 16, 16, 16, 16,
1, 1, 1, 1, 4, 5, 7, 8, 10, 10, 12, 16, 16, 16)
censura=c(0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
1,1,1,0,0,1,1,1,1,0,0,0,0,0)
grupo=c("Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle",
"Controle","Controle",
"Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide",
"Esteróide")

survdiff(Surv(tempos,censura)~grupo,rho=0) #Teste de LogRank (rho=0)
1-pchisq(3.67, df=1, ncp = 0, lower.tail = TRUE, log.p = FALSE)# valor p
qchisq(0.95, df=1, ncp = 0, lower.tail = TRUE, log.p = FALSE) # limiar crítico com alfa=0,05
qchisq(0.9, df=1, ncp = 0, lower.tail = TRUE, log.p = FALSE) # limiar crítico com alfa=0,1


#LOGRANK GENERALIZADO
tempos=c(7,8,8,8,8,12,12,17,18,22,30,30,30,30,30,30,
8,8,9,10,10,14,15,15,18,19,21,22,22,23,25,
8,8,8,8,8,8,9,10,10,10,11,17,19)
censura=c(rep(1,10),rep(0,6),rep(1,15),rep(1,13))
grupo=c(rep("I",16),rep("II",15),rep("III",13))

ekm <- survfit(Surv(tempos,censura)~grupo)
summary(ekm)
autoplot(ekm)
survdiff(Surv(tempos,censura)~grupo,rho=0)
survdiff(Surv(tempos[1:31],censura[1:31])~grupo[1:31],rho=0)#IxII
survdiff(Surv(c(tempos[1:16],tempos[32:44]),c(censura[1:16],censura[32:44]))~c(grupo[1:16],grupo[32:44]),rho=0)#IIxIII
survdiff(Surv(tempos[17:44],censura[17:44])~grupo[17:44],rho=0)#IIxIII

#WILCOXON, TARONE E WARE
packages <- c("survival","ggplot2","ggfortify","survMisc")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
require("survival")
require("ggplot2")
require("ggfortify")
require("survMisc")
require("coin")

tempos=c(1, 2, 3, 3, 3, 5, 5, 16,16,16, 16, 16, 16, 16, 16,
1, 1, 1, 1, 4, 5, 7, 8, 10, 10, 12, 16, 16, 16)
censura=c(0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
1,1,1,0,0,1,1,1,1,0,0,0,0,0)
grupo=c("Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle","Controle",
"Controle","Controle",
"Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide","Esteróide",
"Esteróide")
survdiff(Surv(tempos,censura)~grupo,rho=0)#LOGRANK
survdiff(Surv(tempos,censura)~grupo,rho=1)#WILCOXON
survdiff(Surv(tempos,censura)~grupo,rho=0.5)#TARONE-WARE

# Mais detalhes sobre os testes podem ser encontrados em
# http://www.math.ucsd.edu/~rxu/math284/slect3.pdf

#outras possibilidades
dados=as.data.frame(cbind(tempos,censura,grupo))
dados$tempos=as.numeric(dados$tempos)
dados$censura=as.numeric(dados$censura)
fit = survfit(Surv(tempos, censura) ~ grupo,data=dados)
comp(ten(fit))
# A biblioteca “coin” possui a função logrank_test, que também faz esses testes
# Outra biblioteca útil é a seguinte
#install.packages("survminer")
require("survminer")
ggsurvplot(fit,data=dados, log.rank.weights = "1", pval = TRUE, pval.method = TRUE)
ggsurvplot(fit,data=dados, log.rank.weights = "n", pval = TRUE, pval.method = TRUE)
ggsurvplot(fit,data=dados, log.rank.weights = "sqrtN", pval = TRUE, pval.method = TRUE)
#https://cran.r-project.org/web/packages/survminer/vignettes/Specifiying_weights_in_log-rank_comparisons.html#tharone-ware

