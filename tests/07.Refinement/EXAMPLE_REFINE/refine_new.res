refine
reset
#
#
data xy, noise_osci.dat
#
@newpara_new.mac
#
set cycle,         20
set conver, status:on, dchi:     0.10000000E-001, conf:     0.95000000, chisq:      1.5000000
set relax, start:      1.0000000
#
@k_calc.mac
branch kuplot
  @k_inter.mac
#
run k_calc.mac, plot:k_inter.mac
#
exit
