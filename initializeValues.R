library(data.table)
library(plyr)
library(readxl)

casesCat.dt <- fread("data/catalunya_diari_total_pob.csv",header=TRUE)

casesCat.dt <- casesCat.dt[casesCat.dt$SEXE != "Altres",]
casesCat.dt <- casesCat.dt[casesCat.dt$GRUP_EDAT != "80 o m\xe9s",]
ageGroups <- unique(casesCat.dt$GRUP_EDAT)

startDate = as.Date("2021-05-01")

getCases <- function(n_vacc, n_pre){
  casesCatInter.dt <- fread("data/catalunya_diari_total_pob.csv",header=TRUE)
  casesCatInter.dt <- casesCatInter.dt[casesCatInter.dt$SEXE != "Altres",]
  casesCatInter.dt <- casesCatInter.dt[casesCatInter.dt$GRUP_EDAT != "80 o m\xe9s",]
  
  casesCat.dt <- as.data.table(casesCatInter.dt  |>
                                 dplyr::group_by(GRUP_EDAT, DATA) |>
                                 dplyr::summarise(cases = sum(CASOS_CONFIRMAT),
                                                  hosp = sum(INGRESSOS_TOTAL),
                                                  uci = sum(INGRESSOS_CRITIC),
                                                  death = sum(EXITUS)))
  
  casesCatAgg.dt <-  as.data.table(casesCatInter.dt  |>
                                     dplyr::group_by(DATA) |>
                                     dplyr::summarise(cases = sum(CASOS_CONFIRMAT),
                                                      hosp = sum(INGRESSOS_TOTAL),
                                                      uci = sum(INGRESSOS_CRITIC),
                                                      death = sum(EXITUS)))
  
  
  ageGroups <- unique(casesCat.dt$GRUP_EDAT)
  M = length(ageGroups)
  
  startDate = as.Date("2021-05-01")
  startDatePre = as.Date("2021-05-01") - n_pre
  frac = read.table("data/sequenceFraction.csv")
  frac = frac$fracDelta
  
  n_days <- length(casesCatAgg.dt[casesCatAgg.dt$DATA >= startDate,]$cases) - n_vacc
  n_variant <- min(length(frac), n_days)
  frac <- frac[1:n_variant]
  
  fracBig = rep(1.0, n_days)
  fracBig[1:n_variant] = frac
  
  casesD = as.integer(casesCatAgg.dt[casesCatAgg.dt$DATA >= startDate]$cases[1:n_variant]*frac[1:n_variant]/100)
  cases = matrix(0, ncol = n_days, nrow = M)
  cases_pre = matrix(0, ncol = n_days, nrow = M)
  casesA = matrix(0, ncol = n_days, nrow = M)
  casesB = matrix(0, ncol = n_days, nrow = M)
  casesNoVacc = matrix(0, ncol = n_days, nrow = M)
  hosp = matrix(0, ncol = n_days, nrow = M)
  death = matrix(0, ncol = n_days, nrow = M)
  uci = matrix(0, ncol = n_days, nrow = M)
  
  for(i in 1:M){
    cases[i,] = casesCat.dt[casesCat.dt$GRUP_EDAT == ageGroups[i] &
                              casesCat.dt$DATA >= startDate]$cases[1:n_days]
    cases_pre[i,] = casesCat.dt[casesCat.dt$GRUP_EDAT == ageGroups[i] &
                                  casesCat.dt$DATA >= startDatePre]$cases[1:n_days]
    casesA[i,] = as.integer(cases[i,]*(1.0 - fracBig))
    casesB[i,] = as.integer(cases[i,]*fracBig)
    
    hosp[i,] = casesCat.dt[casesCat.dt$GRUP_EDAT == ageGroups[i] &
                             casesCat.dt$DATA >= startDate]$hosp[1:n_days]
    death[i,] = casesCat.dt[casesCat.dt$GRUP_EDAT == ageGroups[i] &
                              casesCat.dt$DATA >= startDate]$death[1:n_days]
    uci[i,] = casesCat.dt[casesCat.dt$GRUP_EDAT == ageGroups[i] &
                            casesCat.dt$DATA >= startDate]$uci[1:n_days]
  }
  
  casesAgg = as.integer(casesCatAgg.dt[casesCatAgg.dt$DATA >= startDate]$cases[1:n_days])
  hospAgg <- as.integer(casesCatAgg.dt[casesCatAgg.dt$DATA >= startDate]$hosp[1:n_days])
  uciAgg <- as.integer(casesCatAgg.dt[casesCatAgg.dt$DATA >= startDate]$uci[1:n_days])
  deathAgg <- as.integer(casesCatAgg.dt[casesCatAgg.dt$DATA >= startDate]$death[1:n_days])
  
  return(list(casesD = casesD, cases = cases, M = M, hosp = hosp,
              n_days = n_days, n_variant = n_variant, fractionD = frac,
              casesAgg = casesAgg, hospAgg = hospAgg, uciAgg = uciAgg, deathAgg = deathAgg,
              death = death, uci = uci,
              casesA = casesA, casesB = casesB,
              cases_pre = cases_pre))
}

getVaccination<- function(n_pre, n_days, n_vacc){

  casesCat.dt <- fread("data/catalunya_diari_total_pob.csv",header=TRUE)
  casesCat.dt <- casesCat.dt[casesCat.dt$SEXE != "Altres",]
  casesCat.dt <- casesCat.dt[casesCat.dt$GRUP_EDAT != "80 o m\xe9s",]
  vacCat.dt <- as.data.table(casesCat.dt  |>
                               dplyr::group_by(GRUP_EDAT, DATA) |>
                               dplyr::summarise(vac1 = sum(VACUNATS_DOSI_1),
                                                vac2 = sum(VACUNATS_DOSI_2)))
  
  
  ageGroups <- unique(casesCat.dt$GRUP_EDAT)
  M = length(ageGroups)
  
  startDate = as.Date("2021-05-01") - n_pre - n_vacc
  
  v1 <- vector(length = M)
  v2 <- vector(length = M)
  nu1 <- matrix(0, nrow = M, ncol = n_days+n_pre)
  nu2 <- matrix(0, nrow = M, ncol = n_days+n_pre)
  
  for(i in 1:M){
    v2[i] <- sum(vacCat.dt[vacCat.dt$GRUP_EDAT == ageGroups[i] & vacCat.dt$DATA < startDate,]$vac2)
    v1[i] <- sum(vacCat.dt[vacCat.dt$GRUP_EDAT == ageGroups[i] & vacCat.dt$DATA < startDate,]$vac1) - v2[i]
    
    nu2[i,] <- c(vacCat.dt[vacCat.dt$GRUP_EDAT == ageGroups[i] & vacCat.dt$DATA >= startDate,]$vac2[1:(n_days+n_pre-n_vacc)],rep(0, n_vacc))
    nu1[i,] <- c(vacCat.dt[vacCat.dt$GRUP_EDAT == ageGroups[i] & vacCat.dt$DATA >= startDate,]$vac1[1:(n_days+n_pre-n_vacc)],rep(0, n_vacc))
  }
  
  return(list(v1 = v1, v2 = v2, nu1 = nu1, nu2 = nu2, M = M))
}

getContactMatrix <- function(M){
  ageS <- c(4.5, 5.2, 5.6, 5.2, 5.1, 5.7, 6.0, 7.1, 8.5, 8.3, 7.4, 6.6, 5.8, 5.0, 4.6, 3.6)/100
  dataC.dt <- read_excel("data/MUestimates_all_locations_2.xlsx", sheet = "Spain", col_names = F)
  dataC <- as.matrix(dataC.dt)

  dataCNew <- matrix(0.0, nrow = M, ncol = M)
  for(i in 1:M){
    for(j in 1:M){
      dataCNew[i,j] = (ageS[2*(i-1)+1]*(dataC[2*(i-1)+1,2*(j-1)+1] + dataC[2*(i-1)+1,2*(j-1)+2]) +
                         ageS[2*(i-1)+2]*(dataC[2*(i-1)+2,2*(j-1)+1] + dataC[2*(i-1)+2,2*(j-1)+2])) /
        (ageS[2*(i-1)+1] + ageS[2*(i-1)+2])
    }
  }
  
  for(i in 1:M){
    dataCNew[i,] <- dataCNew[i,]/sum(dataCNew[i,])
  }
  C <- dataCNew
  
  ageN <- vector(length = M)
  for(i in 1:M){
    ageN[i] = ageS[2*(i-1)+1] + ageS[2*(i-1)+2]
  }
  NTot <-7.566e6
  N <- ageN * NTot
  
  return(list(N = N, C = C, NTot = NTot))
}

getIHR <- function(){
  ageS <- c(4.5, 5.2, 5.6, 5.2, 5.1, 5.7, 6.0, 7.1, 8.5, 8.3, 7.4, 6.6, 5.8, 5.0, 4.6, 3.6)/100
  ageS <- ageS/sum(ageS)
  IHR <- c(3.0, 0.26, 0.084, 0.042, 0.080, 0.26, 0.40, 0.63, 1.2, 1.9, 2.3, 4.0, 9.6, 10.0, 24.0, 30.5)/100
  pICU <- c(0.243, 0.289, 0.338, 0.389, 0.443, 0.503, 0.57, 0.653, 0.756, 0.866, 0.954, 1.00, 0.972, 
            0.854, 0.654, 0.402)
  pDICU <- c(0.282, 0.286,0.291,0.299,0.310,0.328,0.353,0.390,0.446,0.520,0.604,0.705,0.806, 0.899, 0.969, 1.000)
  pDH <- c(0.039, 0.037, 0.035, 0.035, 0.036, 0.039, 0.045, 0.055, 0.074, 0.107, 0.157, 0.238, 0.353, 0.502, 0.675, 0.832)/0.0832
  
  IHRshort <- c(1:(length(IHR)/2))
  pICUshort <- c(1:(length(IHR)/2))
  pDICUShort <- c(1:(length(IHR)/2))
  pDHShort <- c(1:(length(IHR)/2))
  
  for(i in 1:(length(IHR)/2)){
    IHRshort[i] = (IHR[2*(i-1)+1]*ageS[[2*(i-1)+1]] + IHR[2*(i-1)+2]*ageS[[2*(i-1)+2]])/(ageS[[2*(i-1)+1]] + ageS[[2*(i-1)+2]])
    pICUshort[i] = (pICU[2*(i-1)+1]*ageS[[2*(i-1)+1]] + pICU[2*(i-1)+2]*ageS[[2*(i-1)+2]])/(ageS[[2*(i-1)+1]] + ageS[[2*(i-1)+2]])
    pDICUShort[i] = (pDICU[2*(i-1)+1]*ageS[[2*(i-1)+1]] + pDICU[2*(i-1)+2]*ageS[[2*(i-1)+2]])/(ageS[[2*(i-1)+1]] + ageS[[2*(i-1)+2]])
    pDHShort[i] = (pDH[2*(i-1)+1]*ageS[[2*(i-1)+1]] + pDH[2*(i-1)+2]*ageS[[2*(i-1)+2]])/(ageS[[2*(i-1)+1]] + ageS[[2*(i-1)+2]])
  }
  
  return(list(IHR = IHRshort, pICU = pICUshort, pDICU = pDICUShort, pDH = pDHShort))
}

data = getIHR()
IHR = data$IHR
pICU = data$pICU
pDICU = data$pDICU
pDH = data$pDH


n_vacc <- 15
n_pre <- as.integer(startDate - as.Date("2021-04-01"))
data <- getCases(0, n_pre)
M <- data$M
cases <- data$cases
hosp <- data$hosp
death <- data$death
uci <- data$uci
casesD <- data$casesD
casesAgg <- data$casesAgg
hospAgg <- data$hospAgg
uciAgg <- data$uciAgg
deathAgg <- data$deathAgg
n_days <- data$n_days
n_variant <- data$n_variant
fractionD <- data$fractionD
casesA = data$casesA
casesB = data$casesB
casesPre = data$cases_pre

frac = read.table("data/sequenceFraction.csv")
n_seq = frac$total[1:n_days]
delta_seq = frac$n_delta[1:n_days]


data <- getVaccination(n_pre, n_days, n_vacc)
v1 <- data$v1
v2 <- data$v2
nu1 <- data$nu1
nu2 <- data$nu2

data <- getContactMatrix(M)
N <- data$N
C <- data$C

NTot <- data$NTot
t_boostYoung <- as.integer(as.Date("2021-06-21") - startDate + 1) + n_pre
t_boostYoungClose <- as.integer(as.Date("2021-07-09") - startDate + 1) + n_pre

#parameters 
gA1 = 1.0 - 0.45
gA2 = 1.0 - 0.45
gB1 = 1.0 - 0.4
gB2 = 1.0 - 0.4
sB1 = 1.0 - 0.307
sB2 = 1.0 - 0.88
sA1 = 1.0 - 0.487
sA2 = 1.0 - 0.937

oneDoseAlpha <- (0.8 - (1-sA1))/(1-(1-sA1))
oneDoseDelta <- (0.8 - (1-sB1))/(1-(1-sB1))
twoDoseAlpha <- (0.95 - (1-sA2))/(1-(1-sA2))
twoDoseDelta <- (0.95 - (1-sB2))/(1-(1-sB2))


muA = 1.0/(2.6)
sigmaLatentA = 1.0/(2.6)
muB = 1.0/(2.6)
sigmaLatentB = 1.0/(2.6)

sIHRAlpha <- 1.42
sIHRDelta <- 1.85

stepsDif = 5.0

f_init = rep(0.0, M)
v2[1] <- 1


initialDistCases = apply(casesPre[,1:7], 1, mean)

var1 = 0.037
var2 = 0.12

# data for Stan
data_fit<- list(n_days = n_days, n_variant = n_variant, n_pre = n_pre,  N = N, C = C, M = M, stepsDif = stepsDif, NTot = NTot, 
                cases = cases, hosp_cases = hosp, ICUAgg = uciAgg, DeathAgg = deathAgg, 
                nu1 = nu1, nu2 = nu2, v1 = v1, v2 = v2, 
                gA1 = gA1, gA2 = gA2, gB1 = gB1, gB2 = gB2, 
                sA1 = sA1, sA2 = sA2, sB1 = sB1, sB2 = sB2, 
                oneDoseAlpha = oneDoseAlpha, twoDoseAlpha = twoDoseAlpha, oneDoseDelta = oneDoseDelta, twoDoseDelta = twoDoseDelta,
                muA = muA, sigmaLatentA = sigmaLatentA, muB = muB, sigmaLatentB = sigmaLatentB, 
                sIHRAlpha = sIHRAlpha, sIHRDelta = sIHRDelta,
                IHR = IHR, pICU = pICU, pDICU = pDICU, pDH = pDH,
                t_boostYoung = t_boostYoung, t_boostYoungClose = t_boostYoungClose,
                var1 = var1, var2 = var2,
                f_init = f_init, 
                delta_seq = delta_seq, n_seq = n_seq,
                initialDistCases = initialDistCases)
