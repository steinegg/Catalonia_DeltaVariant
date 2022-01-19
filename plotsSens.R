source("initializeValues.R")
library(cowplot)
library(ggplot2)
library(surveillance)
library(EpiEstim)

ages <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
c_posterior = "#2ca25f"
dateSeq <- seq(from = startDate, by = "day", length = n_days)

#---------------------------  Immunity ----------------------------

fit_stan <- readRDS(file = "data/fit_stan_imm_bin.rds")

#-------------------------------------- Panel 1 -----------------------------------------------------
bs = 12
col = "black"
f = "grey"
s = 1.5
stro = 0.8
alph = 0.5

casAAgg <- apply(casesA, 2, sum)
casBAgg <- apply(casesB, 2, sum)

smr_pred <- cbind(as.data.frame(summary(
  fit_stan, pars = " predCasesD", probs = c(0.025, 0.5, 0.975))$summary))
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names
smr_pred$cases <- casAAgg+casBAgg
smr_pred$date <- seq(from = startDate, by = "day", length = length(smr_pred$mean))

p1 <- ggplot(smr_pred, mapping = aes(x = date)) +
  geom_point(mapping = aes(y = cases),
             shape = 21, colour = col, fill = f, size = s, stroke = stro, alpha = alph) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = c_posterior, alpha = 0.35) +
  geom_line(mapping = aes(y = X50.), color = c_posterior) + 
  labs(x = "Day", y = "inc") +
  xlab("") +
  ylab("Incidence") +
  theme_classic(base_size = bs) +
  theme(axis.text.x = element_blank())

smr_pred <- cbind(as.data.frame(summary(
  fit_stan, pars = "fractionDMod", probs = c(0.025, 0.5, 0.975))$summary))
smr_pred <- smr_pred |>
  dplyr::slice((n_pre+1):length(smr_pred$mean))
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names

smr_pred$cases <- c(fractionD, rep(NA, n_days - n_variant))
smr_pred$date <- seq(from = startDate, by = "day", length = length(smr_pred$mean))

confs = fread("data/conf_seq.csv")
smr_pred$casesUpper <- c(confs$`97.5%`[1:n_days], rep(NA, n_days - n_variant))
smr_pred$casesLower <- c(confs$`2.5%`[1:n_days], rep(NA, n_days - n_variant))

p2 <- ggplot(smr_pred, mapping = aes(x = date)) +
  annotate("rect", fill = "grey", alpha = 0.4, 
           xmin = as.Date("2021-06-20"), xmax =  as.Date("2021-07-09"),
           ymin = -Inf, ymax = Inf) +
  geom_errorbar(aes(ymin=casesLower, ymax=casesUpper), width=.01, alpha = 0.3, size = 0.3,
                position=position_dodge(.1)) +
  geom_point(mapping = aes(y = cases),
             shape = 21, colour = col, fill = f, size = s, stroke = stro, alpha = alph) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = c_posterior, alpha = 0.35) +
  geom_line(mapping = aes(y = X50.), color = c_posterior) + 
  scale_y_continuous(labels=scales::percent) +
  labs(x = "Day", y = "Number of students in bed") +
  coord_cartesian() +
  theme_classic(base_size = bs) +
  xlab("") +
  ylab("% Delta variant") +
  theme(axis.text.x = element_blank())
#ggsave("plots/fitDelta.pdf", device = cairo_pdf)

smr_pred <- cbind(as.data.frame(summary(
  fit_stan, pars = "hospPred", probs = c(0.025, 0.5, 0.975))$summary))
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names
smr_pred$cases <- hospAgg
smr_pred$date <- seq(from = startDate, by = "day", length = length(hospAgg))

p3 <- ggplot(smr_pred, mapping = aes(x = date)) +
  geom_point(mapping = aes(y = cases),
             shape = 21, colour = col, fill = f, size = s, stroke = stro, alpha = alph) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = c_posterior, alpha = 0.35) +
  geom_line(mapping = aes(y = X50.), color = c_posterior) + 
  labs(x = "Day", y = "inc") +
  xlab("") +
  ylab("Hospitalizations") +
  theme_classic(base_size = bs) +
  theme(axis.text.x = element_blank())
#ggsave("plots/dailyHosp.pdf", device = cairo_pdf)

smr_pred <- cbind(as.data.frame(summary(
  fit_stan, pars = "ICUTot", probs = c(0.025, 0.5, 0.975))$summary))
smr_pred <- smr_pred[(n_pre+1):length(smr_pred$mean),]
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names
smr_pred$cases <- uciAgg
smr_pred$date <- seq(from = startDate, by = "day", length = length(hospAgg))

p4 <- ggplot(smr_pred, mapping = aes(x = date)) +
  geom_point(mapping = aes(y = cases),
             shape = 21, colour = col, fill = f, size = s, stroke = stro, alpha = alph) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = c_posterior, alpha = 0.35) +
  geom_line(mapping = aes(y = X50.), color = c_posterior) + 
  #coord_cartesian(ylim = c(0,1000)) +
  labs(x = "Day", y = "inc") +
  xlab("") +
  ylab("ICU admissions") +
  theme_classic(base_size = bs)
#ggsave("plots/dailyICU.pdf", device = cairo_pdf)

smr_pred <- cbind(as.data.frame(summary(
  fit_stan, pars = "DeathTot", probs = c(0.025, 0.5, 0.975))$summary))
smr_pred <- smr_pred[(n_pre+1):length(smr_pred$mean),]
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names
smr_pred$cases <- deathAgg
smr_pred$date <- seq(from = startDate, by = "day", length = length(hospAgg))

p5 <- ggplot(smr_pred, mapping = aes(x = date)) +
  geom_point(mapping = aes(y = cases),
             shape = 21, colour = col, fill = f, size = s, stroke = stro, alpha = alph) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = c_posterior, alpha = 0.35) +
  geom_line(mapping = aes(y = X50.), color = c_posterior) + 
  #coord_cartesian(ylim = c(0,1000)) +
  labs(x = "Day", y = "inc") +
  xlab("") +
  ylab("Fatalities") +
  theme_classic(base_size = bs)
#ggsave("plots/dailyDeaths.pdf", device = cairo_pdf)

smr_pred <- cbind(as.data.frame(summary(
  fit_stan, pars = "incTotNoVaccMod", probs = c(0.025, 0.5, 0.975))$summary))
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names
smr_pred <- smr_pred[(n_pre+1):length(smr_pred$mean),]
smr_pred$date <- seq(from = startDate, by = "day", length = length(smr_pred$X97.5.))

p6 <- ggplot(smr_pred, mapping = aes(x = date)) +
  annotate("rect", fill = "grey", alpha = 0.4, 
           xmin = as.Date("2021-06-20"), xmax =  as.Date("2021-07-09"),
           ymin = -Inf, ymax = Inf) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.,), fill = c_posterior, alpha = 0.35) +
  geom_line(mapping = aes(y = X50.), color = c_posterior) +
  scale_y_continuous(labels=scales::percent) +
  #geom_point(mapping = aes(y = cases),
  #shape = 21, colour = "black", fill = "white", size = 1.3, stroke = 1) +
  #coord_cartesian(ylim = c(0,1000)) +
  labs(x = "Day", y = "inc") +
  xlab("") +
  ylab("% inf. not vaccinated") +
  theme_classic(base_size = bs)

rel_height = 0.93
p <- plot_grid(p1, p3, p2, p4, p5, p6, labels = "AUTO", align = "v", label_size = 12,
               rel_heights = c(rel_height, 1.0))

save_plot("plots/panel1_immun.pdf", p,
          base_height = 5.0,
          base_width = 8.25)

#--------------------------------------------- Panel 2 -----------------------------------------------
bs = 12

smr_pred1 <- cbind(as.data.frame(summary(
  fit_stan, pars = "beta", probs = c(0.025, 0.5, 0.975))$summary))
colnames(smr_pred1) <- make.names(colnames(smr_pred1))
smr_pred1$ages <- rep(ages, each = n_days)
dateVector <- seq(from = startDate, by = "day", length = M*(n_days))
dateSeq <- seq(from = startDate, by = "day", length = n_days)
casesVector = hosp[1,]
index <- 1
for(j in 2:M){
  casesVector <- c(casesVector, hosp[j,])
}
smr_pred1$cases <- casesVector
smr_pred1$date <- rep(dateSeq, M)

smr_pred <- cbind(as.data.frame(summary(
  fit_stan, pars = "repNumber", probs = c(0.025, 0.5, 0.975))$summary))
colnames(smr_pred) <- make.names(colnames(smr_pred))
smr_pred$ages <- rep(ages, each = n_days)
dateVector <- seq(from = startDate, by = "day", length = M*n_days)
dateSeq <- seq(from = startDate, by = "day", length = n_days)
casesVector = hosp[1,]
for(j in 2:M){
  casesVector <- c(casesVector, hosp[j,])
}
smr_pred$cases <- casesVector
smr_pred$date <- rep(dateSeq, M)

p1 <- ggplot(smr_pred1[smr_pred$ages %in% ages[2:6],], mapping = aes(x = date)) +
  geom_vline(xintercept = as.Date("2021-05-09"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-07-15"), linetype = "dashed") +
  annotate("rect", fill = "grey", alpha = 0.4, 
           xmin = as.Date("2021-06-20"), xmax =  as.Date("2021-07-09"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.7, 
           xmin = as.Date("2021-06-23"), xmax =  as.Date("2021-06-27"),
           ymin = -Inf, ymax = Inf) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = c_posterior, alpha = 0.35) +
  geom_line(mapping = aes(y = X50.), color = c_posterior) + 
  facet_wrap(~ ages, ncol = 8) +
  xlab(" ") +
  ylab("Interaction rate") +
  theme_minimal(base_size = bs) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "black", 
                                   size = 0.8, linetype = "solid"),
        axis.line.y = element_line(colour = "black", 
                                   size = 0.8, linetype = "solid"),
        strip.background = element_blank(),
        axis.ticks = element_line(size = 0.5, color="black"),
        axis.text.x = element_blank())

p2 <- ggplot(smr_pred[smr_pred$ages %in% ages[2:6],], mapping = aes(x = date)) +
  geom_vline(xintercept = as.Date("2021-05-09"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-07-15"), linetype = "dashed") +
  annotate("rect", fill = "grey", alpha = 0.4, 
           xmin = as.Date("2021-06-20"), xmax =  as.Date("2021-07-09"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.7, 
           xmin = as.Date("2021-06-23"), xmax =  as.Date("2021-06-27"),
           ymin = -Inf, ymax = Inf) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "#e34a33", alpha = 0.35) +
  geom_line(mapping = aes(y = X50.), color = "#e34a33") + 
  facet_wrap(~ ages, ncol = 8) +
  xlab(" ") +
  ylab("Rep. number") +
  theme_minimal(base_size = bs) +
  scale_y_continuous() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "black", 
                                   size = 0.8, linetype = "solid"),
        axis.line.y = element_line(colour = "black", 
                                   size = 0.8, linetype = "solid"),
        strip.background = element_blank(),
        axis.ticks = element_line(size = 0.5, color="black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        strip.text.x = element_blank())

p <- plot_grid(p1,p2, nrow = 2, align = "v", labels = "AUTO", label_size = 12)
save_plot("plots/panel2_immun.pdf", p,
          base_height = 4.5,
          base_width = 8.25)


#----------------------------------- Panel 3 -------------------------------------------------

#ggsave("plots/casesAges.pdf", device = cairo_pdf)

smr_pred <- smr_pred |>
  dplyr::group_by(date) |>  
  dplyr::summarise(casesAgg = sum(X50.))
casesAggModel <- smr_pred$casesAgg

smr_pred <- cbind(as.data.frame(summary(
  fit_stan, pars = " repTotal", probs = c(0.025, 0.5, 0.975))$summary))
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names
smr_pred$date <- seq(from = startDate, by = "day", length = length(smr_pred$X50.))

l = 2.0/2.6
b = 1.0/2.6
a = l - b

ts = seq(from = 0.0, to = 25, by = 1.0)

ft = l*l*b*exp(-b*ts)*(1/a/a*(1.0 - exp(-a*ts))-ts/a*exp(-a*ts))
ft = ft/sum(ft)

casesCatInter.dt <- fread("catalunya_diari_total_pob.csv",header=TRUE)
casesCatInter.dt <- casesCatInter.dt[casesCatInter.dt$SEXE != "Altres",]
casesCatInter.dt <- casesCatInter.dt[casesCatInter.dt$GRUP_EDAT != "80 o m\xe9s",]
casesCatAgg.dt <-  as.data.table(casesCatInter.dt  |>
                                   dplyr::group_by(DATA) |>
                                   dplyr::summarise(cases = sum(CASOS_CONFIRMAT)))

n_init = 30
casesAgg = casesCatAgg.dt[casesCatAgg.dt$DATA >= (as.Date("2021-05-01") - n_init),]$cases[1:(n_init+n_days)]

incidence.sts <- sts(casesAgg)
bpnp.control <- list(k=0,eps=rep(0.005,2),iter.max=rep(250,2),B=-1,hookFun=NULL,verbose=FALSE)
bpnp.control2 <- modifyList(bpnp.control, list(hookFun=NULL,k=6,B=100,eq3a.method="C"))

sts.bp2 <- backprojNP(incidence.sts, incu.pmf=ft, control=bpnp.control2)

casesBackproj <- upperbound(sts.bp2)[,][(n_init+1):n_days]

t_start = 7
t_window = 7
number_of_days = length(casesBackproj)

config <- make_config(list(si_distr = ft,
                           t_start = seq(t_start, number_of_days+1-t_window),
                           t_end = seq(t_start+t_window-1, number_of_days)))

res.fit.rt <- estimate_R(casesBackproj,
                         method="non_parametric_si",
                         config = config)

res.fit.rtModel <- estimate_R(casesAggModel,
                              method="non_parametric_si",
                              config = config)


Dt = data.table(
  Rt = res.fit.rt$R[3][,1],
  Rmin = res.fit.rt$R[5][,1],
  Rmax = res.fit.rt$R[10][,1],
  RtModel = res.fit.rtModel$R[3][,1],
  RminModel = res.fit.rtModel$R[5][,1],
  RmaxModel = res.fit.rtModel$R[10][,1],
  dates = seq(startDate+(t_window-1)/2+t_start-1, by = "day", length.out = length(res.fit.rt$R[3][,1])))

t1 <- as.Date("2021-06-21") 
t2 <- as.Date("2021-07-09")
Dt1 <- Dt[Dt$dates >= t1 & Dt$dates <= t2,]
Dt2 <- Dt[Dt$dates < t1 || Dt$dates > t2,]

colors <- c("EpiEstim (Data)" = "black", "Model" = "red")

p1 <- ggplot() +
  geom_vline(xintercept = as.Date("2021-05-09"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-07-15"), linetype = "dashed") +
  annotate("rect", fill = "grey", alpha = 0.4, 
           xmin = as.Date("2021-06-20"), xmax =  as.Date("2021-07-09"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.7, 
           xmin = as.Date("2021-06-23"), xmax =  as.Date("2021-06-27"),
           ymin = -Inf, ymax = Inf) +
  scale_color_manual(values = colors) +
  geom_ribbon(data = Dt, aes(x = dates, ymin = Rmin, ymax = Rmax), color = "black", alpha = 0.35) +
  geom_line(data = Dt, mapping = aes(x = dates, y = Rt, color = "EpiEstim (Data)")) +
  geom_point(data = Dt, mapping = aes(x = dates, y = Rt), color = "black", shape = 18, size = 2) +
  ylab("Rep. number") +
  xlab(" ") +
  labs(color = "") +
  geom_ribbon(data = smr_pred, aes(x= date, ymin = X2.5., ymax = X97.5.), fill = "red", alpha = 0.35) +
  geom_line(data = smr_pred, mapping = aes(x = date, y = X50.), color = "red") +
  theme_classic(base_size = bs) +
  theme(legend.position = c(0.75, 0.5))

smr_pred <- cbind(as.data.frame(summary(
  fit_stan, pars = " fDTime", probs = c(0.025, 0.5, 0.975))$summary))
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names
smr_pred$date <- seq(from = startDate, by = "day", length = length(smr_pred$mean))

p2 <- ggplot(smr_pred, mapping = aes(x = date)) +
  geom_vline(xintercept = as.Date("2021-05-09"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-07-15"), linetype = "dashed") +
  annotate("rect", fill = "grey", alpha = 0.4, 
           xmin = as.Date("2021-06-20"), xmax =  as.Date("2021-07-09"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.7, 
           xmin = as.Date("2021-06-23"), xmax =  as.Date("2021-06-27"),
           ymin = -Inf, ymax = Inf) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = c_posterior, alpha = 0.35) +
  geom_line(mapping = aes(y = X50.), color = c_posterior) + 
  labs(x = "Day", y = "inc") +
  xlab("") +
  ylab("Detection rate") +
  theme_classic(base_size = bs)
#ggsave("plots/dailyinc.pdf", device = cairo_pdf)

smr_pred <- cbind(as.data.frame(summary(
  fit_stan, pars = "attackRate", probs = c(0.025, 0.5, 0.975))$summary))
colnames(smr_pred) <- make.names(colnames(smr_pred))
smr_pred$ages <- ages
attackRates <- apply(cases[,], 1, sum)/N
smr_pred$reported <- attackRates

p4 <- ggplot(smr_pred) +
  geom_errorbarh(aes(xmin=X2.5., xmax=X97.5., y = ages), 
                 position = position_dodge(0.3), height=0.0, color=c_posterior, 
                 alpha=0.35, size=2) +
  geom_point(aes(x = X50., y = ages), shape = 5, size = 1.5, color = c_posterior, stroke = 1.1) +
  scale_x_continuous(labels = scales::percent, breaks = c(0.0, 0.15, 0.3, 0.45)) +
  ylab(" ") +
  xlab("Attack rate model") +
  theme_classic(base_size = bs) +
  theme(axis.text.y = element_blank())

p5 <- ggplot(smr_pred) +
  geom_point(aes(x = attackRates, y = ages), shape = 5, size = 1.5, color = "#2b8cbe", stroke = 1.3) +
  scale_x_continuous(labels = scales::percent) +
  ylab(" ") +
  xlab("Attack rate reported") +
  theme_classic(base_size = bs) +
  theme(axis.text.y = element_blank())

smr_pred <- cbind(as.data.frame(summary(
  fit_stan, pars = "fD", probs = c(0.025, 0.5, 0.975))$summary))
colnames(smr_pred) <- make.names(colnames(smr_pred))
smr_pred$ages <- ages

p3 <- ggplot(smr_pred) +
  geom_errorbarh(aes(xmin=X2.5., xmax=X97.5., y = ages), 
                 position = position_dodge(0.3), height=0.0, color="black", 
                 alpha=0.35, size=2) +
  geom_point(aes(x = X50., y = ages), shape = 5, size = 1.5, color = "black", stroke = 1.1) +
  scale_x_continuous(labels = scales::percent) +
  ylab(" ") +
  xlab("Detection rate") +
  theme_classic(base_size = bs)

p6 <- plot_grid(p1, p2, ncol = 2, labels = c("A", "B"), label_size = 12, align = "h")
p7<- plot_grid(p3, p5, p4, ncol = 3, rel_widths = c(1., .9, .9), labels = c("C", "D", "E"), label_size = 12)

p <- plot_grid(p6, p7, ncol = 1)

save_plot("plots/panel3_immun.pdf", p,
          base_height = 5.0,
          base_width = 8.25)

#--------------------------- Severity ----------------------------
fit_stan <- readRDS(file = "data/fit_stan_sev_bin.rds")

smr_pred <- cbind(as.data.frame(summary(
  fit_stan, pars = "attackRate", probs = c(0.025, 0.5, 0.975))$summary))
colnames(smr_pred) <- make.names(colnames(smr_pred))
smr_pred$ages <- ages
attackRates <- apply(cases[,], 1, sum)/N
smr_pred$reported <- attackRates

p <- ggplot(smr_pred) +
  geom_errorbarh(aes(xmin=X2.5., xmax=X97.5., y = ages), 
                 position = position_dodge(0.3), height=0.0, color=c_posterior, 
                 alpha=0.35, size=2) +
  geom_point(aes(x = X50., y = ages), shape = 5, size = 1.5, color = c_posterior, stroke = 1.1) +
  scale_x_continuous(labels = scales::percent, breaks = c(0.0, 0.15, 0.3, 0.45)) +
  ylab(" ") +
  xlab("Attack rate model") +
  theme_classic(base_size = bs) +
  theme()
save_plot("plots/severity.pdf", p)

#--------------------------- Generation Time ----------------------------
pars=c("fB")
fit_stan <- readRDS(file = "data/fit_stan_Main_bin.rds")
print(fit_stan, pars = pars)

res.dt <- data.table(
  max = c(1.43, 1.54, 1.66, 1.48),
  min = c(1.39, 1.50, 1.60, 1.44),
  med = c(1.41, 1.52, 1.63, 1.46),
  type = c("4.2", "5.2 (main)", "6.2", "4.6/5.5")
)

p <- ggplot(res.dt) +
  geom_linerange(aes(ymin=min, ymax=max, x = type), 
                 position = position_dodge(0.3),  color=c_posterior, 
                 alpha=0.35, size=2) +
  geom_point(aes(x = type, y = med), shape = 5, size = 1.5, color = c_posterior, stroke = 1.1) +
  xlab(" ") +
  ylab("Transm. advantage") +
  theme_classic(base_size = bs) +
  theme()
save_plot("plots/genTime.pdf", p)
