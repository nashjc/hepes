# Observed data: 86 instances with variable SID and HEPES and measured pH from
# the sulfate paper plus 21 with variable SID  fixed HEPES and measured pH.
# Fitted pH obtained for monovalent HEPES with pKa 7.55 using charge-balance modeling
# and with trivalent HEPES with pK = c(-1,3,7.55) using the Kildeberg algorithm and
# (identical) ad hoc formula derived here. Also the formulation of Johnson 2014 with
# pK1 = 3 and pK2 = 7.55
# 
#  First get observed data
load('HEPES.RData') # DD: 86 lines data.frame SID, HEPES, pH obs
BUFC = 0.050
pHobs  <- c(4.63,4.68,4.72,4.77,4.83,4.9,4.96,5.04,5.12,5.21,5.3,
            5.39,5.48,5.55,5.63,5.69,5.74,5.8,5.85,5.89,5.93)

SIDobs<- c(-seq(10,1)*1e-4,0,seq(1,10)*1e-4)
SID <- c(DD$SID,SIDobs)
pH <- c(DD$OBS,pHobs)
HEPES <- c(DD$HEPES,rep(BUFC,21))



kw <- 1e-14


library(nlsr)
library(ggplot2)
library(ggpubr)

HEPESFUNC <- function(H,SID,HEPTOT,pK1,pK2,pK3) {
  XX <- (H^3/(10^-pK3*10^-pK2*10^-pK1)+H^2/(10^-pK3*10^-pK2)+H/(10^-pK3))
  IV <- HEPTOT/(1+XX)
  I <- IV*H^3/(10^-pK3*10^-pK2*10^-pK1)
  II <- IV*H^2/(10^-pK3*10^-pK2)
  III <- IV*H/10^-pK3
  H - kw/H + SID + I*2 + II - IV
}



pHm <- Vectorize(function(SID, HEPTOT, pK1, pK2, pK3)
  -log10(uniroot(HEPESFUNC,c(1e-20,1),tol=1e-20,maxiter=1e4,
                 HEPTOT=HEPTOT,SID = SID, pK1=pK1,pK2=pK2,pK3=pK3)$root))

coefpp <- function(model) {
  co <- coef(model)
  if (identical(names(co), paste0('pK', 1:3))) return(co)
  if (identical(names(co), c('pK1', paste0('d', c('12', '23')))))
    return(setNames(cumsum(co), paste0('pK', 1:3)))
  stop('What do I do with parameters named ', deparse(names(co)), '?')
}


# Kildeberg formula for mean charge on HEPES given ordered vector of pKa;
# n denotes charge on completely acidified species
NN <- function(H,pK,n) {
  L <- 1:length(pK)
  num <- 0
  for (i in 1:length(pK)) num <- num+L[i]*10^(L[i]*(-log10(H))-sum(pK[1:i]))
  denum <- 1
  for (i in 1:length(pK)) denum <- denum+10^(L[i]*(-log10(H))-sum(pK[1:i]))
  num/denum - n
}

# Multiprotic am Kildeberg
BUFFH <- function(H,HEPTOT,SID,pK,n) {
  SID+H-kw/H - HEPTOT*NN(H,pK,n)}

# Monoprotic
HEPmono <- function(SID,H,HEPTOT,pK) { 
  KA <- 10^-pK
  SID+H-kw/H - HEPTOT*KA/(KA+H)}

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
pKg <- c(-1,3,7.55)
plot(SID,pH)

mn <- mu <- c()
for (i in 1:length(SID)) {
  mn[i] <- -log10(uniroot(HEPmono,c(1e-20,1),
                          tol=1e-20,maxiter=1e4,SID=SID[i],HEPTOT=HEPES[i],pK=7.55)$root)
  mu[i] <-  -log10(uniroot(BUFFH,c(1e-20,1),
                           tol=1e-20,maxiter=1e4,SID=SID[i],HEPTOT=HEPES[i],pK=pKg,n=2)$root)
}

points(SID,mn,col="blue")
points(SID,mu,col="red")


(ssq <- sum((pH-mn)^2)) # 26.30414 monoprotic SSQ
(ssq <- sum((pH-mu)^2)) # 0.1381444 multiprotic SSQ
plot(pH,mn)
abline(0,1,col="red") # acid deviation
plot(pH,mu)
abline(0,1,col="red") # more reasonable
ddf <- data.frame(SID=SID,pH=pH,HEPES=HEPES,mn=mn,mu=mu)





#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Use pK 3, 7.5 a.m. Johnson 2014
#  First Kildeberg with pKs 3 and 7.5 and n = 1 
#  corresponding to sulfonate group neutral and
#  piperazine nitrogen +1


pHex <- c()
for (i in 1:length(SID)) {
  pHex[i] <- -log10(uniroot(BUFFH,c(1e-20,1),tol=1e-20,maxiter=1e4,
                            HEPTOT=HEPES[i],SID=SID[i],pK=c(3,7.5),n=1)$root)
}

plot(SID,pH)
points(SID,mu,col="red")
points(SID,pHex,col="blue")

plot(pHex,pH)
abline(0,1,col="blue")

(SSQ <- sum((pHex-pH)^2)) # 0.2299074

# So it is evident that  this specification is inferior - this was using Kildeberg
# Below we use direct evaluation as explained in 
# HEPES_n1_061123.docx

# We have H + SID - kw/H + I - III = 0
# III = HEPTOT/(1+ H^2/(10^-7.5 * 10^-3)+ H/10^-7.5)
# I = III*H^2/(10^-7.5 * 10^-3)

HEPJohns <- function(H,SID,HEPTOT) {H + SID - kw/H - 
    HEPTOT/(1+ H^2/(10^-7.5 * 10^-3)+ H/10^-7.5)+
    HEPTOT/(1 + H^2/(10^-7.5 * 10^-3)+ H/10^-7.5)*H^2/(10^-7.5 * 10^-3)}

pHJohns <- c()
for (i in 1:length(SID)) {
  pHJohns[i] <- -log10(uniroot(HEPJohns,c(1e-20,1),tol=1e-20,maxiter=1e4,
                               HEPTOT=HEPES[i],SID=SID[i])$root)
}
plot(pHJohns,pHex)
abline(0,1,col="red")
all.equal(pHJohns,pHex)
#TRUE

cbind(pHJohns,pHex)
#QED!!!!!!!

ddf <- data.frame(SID=SID,pH=pH,HEPES=HEPES,mn=mn,mu=mu,pHjohns = pHex)



# make a plot plus regression on obs pH vs fitted with conventional pks
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# sum of square  devs for each BA plots
(ssq <- sum((ddf$pH-ddf$mu)^2)) #  0.13814
GG <- ggplot()
GG <- GG + geom_point(data=ddf,aes(x=pH,y=mu),size=2)+
  geom_abline(intercept=0,slope=1,col="red",linewidth=1)+
  theme(axis.title=element_text(size=15))+
  theme(axis.text=element_text(size=15))+
  xlab("Observed pH") + ylab("Modeled pH")+
  ggtitle(label="A: Modeled versus measured pH", 
          subtitle =expression(paste(pK[1]," = -1, ",pK[2]," = 3, ",pK[3]," = 7.55 Red line marks X = Y")))

GG
cor(pH,mu) # 0.99935
summary(z1 <- lm(mu~pH))

# Monoprotic 
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# 
FF <- ggplot()
FF <- FF + geom_point(data=ddf,aes(x=pH,y=mn),size=2)+
  geom_abline(intercept=0,slope=1,col="red",linewidth=1)+
  theme(axis.title=element_text(size=15))+
  theme(axis.text=element_text(size=15))+
  xlab("Observed pH") + ylab("Modeled pH")+
  ggtitle(label="B: Modeled with single pKa versus measured pH", 
          subtitle =expression(paste(pK[1]," = 7.55 Red line marks X = Y")))

FF

HH <- ggplot()
HH <- HH + geom_point(data=ddf,aes(x=pH,y=pHjohns),size=2)+
  geom_abline(intercept=0,slope=1,col="red",linewidth=1)+
  theme(axis.title=element_text(size=15))+
  theme(axis.text=element_text(size=15))+
  xlab("Observed pH") + ylab("Modeled pH")+
  ggtitle(label="C:  Modeled am Johnson versus measured pH", 
          subtitle =expression(paste(pK[1]," = 7.55 Red line marks X = Y")))

HH

(gfh <- ggarrange(GG,FF,HH,nrow=1))
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

saveRDS(gfh, 'observed_final_250124.rds', version = 2, compress = 'xz')

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Fit using nlsr ordering specified in nllsr

HEPTOT <- HEPES

# fixed pK1 = -1 
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
pK1 <- -1; pK2 <- 3; pK3 <- 7.55 # literature values
nlsr::nlxb(
  pHobs ~ pHm(SID, HEPTOT,pK1,pK1+d12, pK1+d12+d23),
  data.frame(pHobs=pH, SID = SID),
  start = c(pK1 = pK1, d12 = pK2 - pK1, d23 = pK3 - pK2),
  # lower bounds provide the ordering
  lower = c(pK1 = -1, d12 = .001, d23 = .001),
  # the original upper bounds on pK may be exceeded in practice
  upper = c(pK1 = -1, d12 = 9 - -3, d23 = 9 - -3),
  trace = TRUE, # so far nlxb hasn't crashed on me yet
  # the following (or some other manual differentiation fix) is needed
  # because nlxb rightfully doesn't know how to differentiate `pHm`
  control = list(japprox = 'jacentral')
) -> m
coefpp(m)
# pK1     pK2     pK3 
# -1.0000  2.9969  7.5350 

ddf$nlsr <- as.numeric(fitted(m))

cor(ddf$pH,ddf$nlsr) # 0.9993359
plot(ddf$pH,ddf$nlsr)
abline(0,1,col="red")
m
# residual sumsquares =  0.117  on  107 observations
# 
(ssq <- sum((ddf$pH-ddf$nlsr)^2)) #  0.117

cor(pH,ddf$nlsr) # 0.99934
summary(z1 <- lm(ddf$nlsr~pH))


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#Ordering ignored in nlsr
#
#pK1 fixed <- -1; pK2 <- 3; pK3 <- 7.55 # literature values
nlsr::nlxb(
  pHobs ~ pHm(SID, HEPTOT,pK1,pK2, pK3),
  data.frame(pHobs=pH, SID = SID),
  start = c(pK1 = pK1, pK2=pK2, pK3=pK3),
  # lower bounds provide the ordering
  lower = c(pK1 = -1, pK2=1, pK3=1),
  # the original upper bounds on pK may be exceeded in practice
  upper = c(pK1 = -1, pK2 = 9 , pK3 = 9),
  trace = TRUE, # so far nlxb hasn't crashed on me yet
  # the following (or some other manual differentiation fix) is needed
  # because nlxb rightfully doesn't know how to differentiate `pHm`
  control = list(japprox = 'jacentral')
) -> m
coefpp(m)
# pK1     pK2     pK3 
# -1.0000  2.9969  7.5350 

ddf$nlsr2 <- as.numeric(fitted(m))
all.equal(ddf$nlsr2,ddf$nlsr)
#  TRUE
