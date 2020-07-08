# cheatsheet for which HH model to plot
# NOTE: this was created done by hand and is utterly shameful
which.model <- as.array(
  rbind(
  c("Chan_2017_c","identity"),
  c("Chan_2017_l","exp"),
  c("Colton_1987_1","identity"),
  c("Colton_1987_2","identity"),
  c("Elliot_2006_i2","identity"),
  c("Elliot_2006_i3","exp"),
  c("Elliot_2006_i4","identity"),
  c("Elliot_2006_i4B","identity"),
  c("Elliot_2006_i5","identity"),
  c("Elliot_2006_i5B","identity"),
  c("Iyer_1996_Bc","identity"),
  c("Iyer_1996_Bp","identity"),
  c("Iyer_1996_Br","identity"),
  c("Kalinkat_2011_A","identity"),
  c("Kalinkat_2011_C","exp"),
  c("Kalinkat_2011_H","identity"),
  c("Kalinkat_2011_P","exp"),
  c("Kalinkat_2011_T","exp"),
  c("Krylov_1992_ii","identity"),
  c("Lester_2002_Af_d","identity"),
  c("Lester_2002_Af_e","identity"),
  c("Lester_2002_Ty_d","identity"),
  c("Lester_2002_Ty_e","identity"),
  c("Long_2012_2p","identity"),
  c("Nachappa_2006","exp"),
  c("Ranta_1985_10","identity"),
  c("Ranta_1985_13","identity"),
  c("Ranta_1985_18","identity"),
  c("Ranta_1985_Ad","exp"),
  c("Wong_2005_rc","identity"),
  c("Wong_2005_ss","exp")
  )
)
rownames(which.model) <- which.model[,1]
which.model <- which.model[,2]