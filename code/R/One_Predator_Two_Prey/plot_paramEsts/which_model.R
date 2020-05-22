# cheatsheet for which HH model to plot
# NOTE: this was created done by hand and is utterly shameful
which.model <- as.array(
  rbind(
  c("Chan_2017_c","identity"),
  c("Chan_2017_l","exp"),
  c("Colton_1987_1st24","identity"),
  c("Colton_1987_2nd24","identity"),
  c("Elliot_2006_Instar2","identity"),
  c("Elliot_2006_Instar3","exp"),
  c("Elliot_2006_Instar4","identity"),
  c("Elliot_2006_Instar4Baet","identity"),
  c("Elliot_2006_Instar5","identity"),
  c("Elliot_2006_Instar5Baet","identity"),
  c("Iyer_1996_Bc","identity"),
  c("Iyer_1996_Bp","identity"),
  c("Iyer_1996_Br","identity"),
  c("Kalinkat_2011_Anch","identity"),
  c("Kalinkat_2011_Cal","exp"),
  c("Kalinkat_2011_Harp","identity"),
  c("Kalinkat_2011_Pard","exp"),
  c("Kalinkat_2011_Troch","exp"),
  c("Krylov_1992ii","identity"),
  c("Lester_2002_Af_duets","identity"),
  c("Lester_2002_Af_eggs","identity"),
  c("Lester_2002_Ty_duets","identity"),
  c("Lester_2002_Ty_eggs","identity"),
  c("Long_2012b","identity"),
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