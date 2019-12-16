# Some simple stats for the Sulman Z* paper
# 2019-11-19

meanSE=function(data,tapplyString1){
  
  data_mean=tapply(data, tapplyString1, mean, na.rm=T)
  data_sd=tapply(data, tapplyString1, sd, na.rm=T)
  data_count=tapply(!is.na(data),tapplyString1, sum)
  data_se= data_sd/sqrt(data_count)
  data_med = tapply(data, tapplyString1, median, na.rm=T)
  list(means=data_mean, se=data_se, sd=data_sd, n=data_count, median = data_med)
}

dat = read.table('profiledata_stats_bins_Ap_20190624_ct.csv', sep=',', header = T)


# Anovas / paired t-tests for the calculated parameters by with and without Ap horizons

# Figure 5. Distribution
# chi2 test for the distribution (2-way table)
table(dat$has_Ap)

# Zstar, 5 cm bins
Fig5 = function(VarIn, datIn, binlim){
  zsT = hist(VarIn[which(datIn$has_Ap == T)], breaks = binlim, freq=T)
  zsF = hist(VarIn[which(datIn$has_Ap == F)], breaks = binlim, freq=T)

  Zs1 = as.table(cbind(zsT$counts, zsF$counts))
  dimnames(Zs1)=list('Zstarval'= zsT$mids, 'ApTF' = c('ApT', 'ApF'))
  Zchisq1 = chisq.test(matrix(Zs1))
  print('Chi squared, Ap vs. none')
  print(Zchisq1)
  
  ZsTm = hist(VarIn[which(datIn$has_Ap == T & datIn$soilorder == 'Mollisol')], breaks = binlim, freq=T)
  ZsFm = hist(VarIn[which(datIn$has_Ap == F & datIn$soilorder == 'Mollisol')], breaks = binlim, freq=T)
  Zchisq2 = chisq.test(matrix(cbind(ZsTm$counts, ZsFm$counts)))
  print('Chi squared, Mollisol: Ap vs none')
  print(Zchisq2)
  
  ZsMnot = hist(VarIn, breaks = binlim, freq=T)
  ZsM = hist(VarIn[which(datIn$soilorder == 'Mollisol')], breaks = binlim, freq=T)
  Zchisq3 = chisq.test(matrix(cbind(ZsMnot$counts, ZsM$counts)))
  print('Chi squared: Mollisol vs all soils')
  print(Zchisq3)

  # find medians & 95% confidence intervals (2*sd) using bootstrap resampling
  nreps =10e3
  Zs_T_med = rep(NA, times=nreps)
  Zs_F_med = rep(NA, times=nreps)
  Zs_Tm_med = rep(NA, times=nreps)
  Zs_Fm_med = rep(NA, times=nreps)
  Zs_All_med = rep(NA, times=nreps)
  Zs_M_med = rep(NA, times=nreps)
  
  for (i in 1:nreps){
    Zs_T_med[i] = median(sample(VarIn[which(datIn$has_Ap==T)], replace = T))
    Zs_F_med[i] = median(sample(VarIn[which(datIn$has_Ap==F)], replace = T))
    Zs_Tm_med[i] = median(sample(VarIn[which(datIn$has_Ap==T & datIn$soilorder == 'Mollisol')], replace = T))
    Zs_Fm_med[i] = median(sample(VarIn[which(datIn$has_Ap==F & datIn$soilorder == 'Mollisol')], replace = T))
    Zs_All_med[i] = median(sample(VarIn, replace = T))
    Zs_M_med[i] = median(sample(VarIn[which(datIn$soilorder == 'Mollisol')], replace = T))
  }
  print(c('Ap Median 95% CI', round(mean(Zs_T_med), 3), round(2*sd(Zs_T_med), 3)))
  print(c('No Ap, Median 95% CI', round(mean(Zs_F_med),2), round(2*sd(Zs_F_med), 3)))
  print(c('Mollisol Ap, median 95% CI', round(mean(Zs_Tm_med),2), round( 2*sd(Zs_Tm_med),3)))
  print(c('Mollison, no Ap, median 95% CI', round(mean(Zs_Fm_med), 3), round(2*sd(Zs_Fm_med), 3)))
  print(c('All soils, median 95% CI', round(mean(Zs_All_med), 3), round(2*sd(Zs_All_med),3)))
  print(c('All mollisols', round(mean(Zs_M_med),3), round(2*sd(Zs_Fm_med),3)))
}

zstar_seq = seq(0, ceiling(max(dat$zstar)/5)*5, by=5)
Fig5(dat$zstar, dat, zstar_seq)
Csurf_seq = seq(0, ceiling(max(dat$Csurf)), by=.01)
Fig5(dat$Csurf, dat, Csurf_seq)
zmin_seq = seq(0, ceiling(max(dat$Zmin)/10)*10, by=10)
Fig5(dat$Zmin, dat, zmin_seq)
Cdeep_seq = seq(0, ceiling(max(dat$Cdeep)/.001)*.001, by=.001)
  Fig5(dat$Cdeep, dat, Cdeep_seq)
  nreps =10e3
  Zs_T_med = rep(NA, times=nreps)
  Zs_F_med = rep(NA, times=nreps)
  for (i in 1:nreps){
    Zs_T_med[i] = median(sample(dat$Cdeep[which(dat$has_Ap==T)], replace = T))
    Zs_F_med[i] = median(sample(dat$Cdeep[which(dat$has_Ap==F)], replace = T))
  }
  print(c('Ap Median 95% CI', mean(Zs_T_med), 2*sd(Zs_T_med)))
  print(c('No Ap, Median 95% CI', mean(Zs_F_med), 2*sd(Zs_F_med)))
  
Cs_seq = seq(0, ceiling(max(dat$Cstocks_zstar_30cm)/.5)*.5, by=0.5)
Fig5(dat$Cstocks_zstar_30cm, dat, Cs_seq)


#### Figure 6
dat$landuse2 = NA
dat$landuse2[which(dat$landuse == 'Cultivated Crops')] <- 'Cultivated Crops'
dat$landuse2[which(dat$landuse == 'Pasture/Hay')]       <- 'Pasture/Hay'
dat$landuse2[which(dat$landuse == 'Grassland/Herbaceous')] <- 'Grassland/Herbaceous'
dat$landuse2[which(dat$landuse == 'Shrub/Scrub')] <- 'Shrub/Scrub'
dat$landuse2[which(dat$landuse == 'Deciduous Forest')] <- 'Deciduous Forest'
dat$landuse2[which(dat$landuse == 'Mixed Forest')] <- 'Mixed Forest'
dat$landuse2[which(dat$landuse == 'Evergreen Forest')] <- 'Evergreen Forest'
dat$landuse2 = as.factor(dat$landuse2)
dat$For_nonFor = as.character(dat$landuse2)
dat$For_nonFor[which(dat$landuse2 == 'Cultivated Crops' | dat$landuse2 == 'Pasture/Hay' |  
                       dat$landuse2 == 'Grassland/Herbaceous' | dat$landuse2 == 'Shrub/Scrub')] <- 'Non-forest'
dat$For_nonFor[which(dat$landuse2 == 'Deciduous Forest' | dat$landuse2 == 'Mixed Forest' |  
                       dat$landuse2 == 'Evergreen Forest')] <- 'Forest'
dat$For_nonFor = as.factor(dat$For_nonFor)

(aov1 = aov(dat$Csurf ~ as.factor(dat$has_Ap) * dat$landuse2))
summary(aov1)
TukeyHSD(aov1)
(aov1A = aov(dat$Csurf ~ as.factor(dat$has_Ap)))
summary(aov1A)
meanSE(dat$Csurf, as.factor(dat$has_Ap))

(aov2 = aov(dat$Cstocks_zstar_30cm ~ as.factor(dat$has_Ap) * dat$landuse2))
summary(aov2)
meanSE(dat$Cstocks_zstar_30cm, as.factor(dat$has_Ap))

(aov3 = aov(dat$Csurf ~ as.factor(dat$has_Ap) * dat$For_nonFor))
summary(aov3)
TukeyHSD(aov3)
meanSE(dat$Csurf, list(as.factor(dat$has_Ap), dat$For_nonFor))

(aov3A = aov(dat$zstar ~ as.factor(dat$has_Ap) * dat$For_nonFor))
summary(aov3A)
TukeyHSD(aov3A)
meanSE(dat$zstar, list(as.factor(dat$has_Ap), dat$For_nonFor))
(aov3B = aov(dat$zstar ~ as.factor(dat$For_nonFor)))
summary(aov3B)


(aov1A = aov(dat$Csurf ~ as.factor(dat$has_Ap)))
summary(aov1A)
meanSE(dat$Csurf, as.factor(dat$has_Ap))

summary(aov(Csurf ~ landuse2, data=dat[which(dat$has_Ap == T), ]))
meanSE(dat$Csurf[which(dat$has_Ap == T)], dat$landuse2[which(dat$has_Ap ==T)])


# Reclass the ones without enough data
dat$soilorderF6 = as.character(dat$soilorder)
dat$soilorderF6[which(dat$soilorder == "Andisol")] <- NA
dat$soilorderF6[which(dat$soilorder == "Gelisol")] <- NA
dat$soilorderF6[which(dat$soilorder == "Histosol")] <- NA
dat$soilorderF6[which(dat$soilorder == "Oxisol")] <- NA
dat$soilorderF6 = as.factor(dat$soilorderF6)

(aov10 = aov(dat$Csurf ~ as.factor(dat$has_Ap) * dat$soilorderF6))
summary(aov10)
TukeyHSD(aov10)
a10 = meanSE(dat$Csurf, list(as.factor(dat$has_Ap), dat$soilorderF6))
list2ascii(a10, 'Csurf_means')
summary((a10$means[2,] - a10$means[1,])/a10$means[1,])

(aov20 = aov(dat$Cstocks_zstar_30cm ~ as.factor(dat$has_Ap) * dat$soilorderF6))
summary(aov20)
TukeyHSD(aov20)
a20 = meanSE(dat$Cstocks_zstar_30cm, list(as.factor(dat$has_Ap), dat$soilorderF6))
a20d = (a20$means[2,] - a20$means[1,])/a20$means[1,]
summary(a20d[-c(2,3,6,9)])


(aov30 = aov(dat$zstar ~ as.factor(dat$has_Ap) * dat$soilorderF6))
summary(aov30)
TukeyHSD(aov30)
a30 = meanSE(dat$zstar, list(as.factor(dat$has_Ap), dat$soilorderF6))
summary((a30$means[2,] - a30$means[1,])/a30$means[1,])




############
# Figure 7

t1 = table(as.factor(dat$cell_number), dat$soilorder)
t1 = tapply(dat$Cstocks_zstar_30cm, list(as.factor(dat$cell_number), dat$soilorder, dat$has_Ap), mean)
t.test(t1[,,1], t1[,,2], paired = T)
t1_noAP = data.frame(t1[,,1])
t1_noAP$CellNo = row.names(t1_noAP)
t1_noAP$CellNo = as.factor(t1_noAP$CellNo)
t1_noAP_l = reshape(t1_noAP, direction = 'long', idvar = 'Cellno', varying = levels(dat$soilorder), v.names = 'Cstock')
t1_AP = data.frame(t1[,,2])
t1_AP$CellNo = row.names(t1_AP)
t1_AP$CellNo = as.factor(t1_AP$CellNo)
t1_AP_l = reshape(t1_AP, direction = 'long', idvar = 'Cellno', varying = levels(dat$soilorder), v.names = 'Cstock')

f7 = t.test(t1_AP_l$Cstock, t1_noAP_l$Cstock, paired=T)
t1_diff = t1_noAP_l$Cstock - t1_AP_l$Cstock
br_t1= seq(floor(min(t1_diff, na.rm=T)), ceiling(max(t1_diff, na.rm=T)), 1)
t1_hist = hist(t1_diff, breaks=br_t1)
