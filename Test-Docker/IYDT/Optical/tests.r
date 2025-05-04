########################################################################
#                                                                      #
# © 2023. Triad National Security, LLC. All rights reserved.           #
# This program was produced under U.S. Government contract             #
# 89233218CNA000001 for Los Alamos National Laboratory (LANL), which   #
# is operated by Triad National Security, LLC for the U.S. Department  #
# of Energy/National Nuclear Security Administration. All rights in    #
# the program are reserved by Triad National Security, LLC, and the    #
# U.S. Department of Energy/National Nuclear Security Administration.  #
# The Government is granted for itself and others acting on its behalf #
# a nonexclusive, paid-up, irrevocable worldwide license in this       #
# material to reproduce, prepare derivative works, distribute copies   #
# to the public, perform publicly and display publicly, and to permit  #
# others to do so.                                                     #
#                                                                      #
########################################################################

# directory to forward models and jacobians
adir = "Code/IYDT-gsrp"

# source forward models and jacobians
source(paste(adir,"/forward.r",sep=""),local=TRUE)
source(paste(adir,"/forward_0.r",sep=""),local=TRUE)
source(paste(adir,"/jacobian.r",sep=""),local=TRUE)
source(paste(adir,"/jacobian_0.r",sep=""),local=TRUE)

# directory to test files
tdir = "Test"

print("***** test_rom_o1.r *****")
cat("\n")
source(paste(tdir,"/test_rom_o1.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_o2-0.r *****")
cat("\n")
source(paste(tdir,"/test_rom_o2-0.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_o2.r *****")
cat("\n")
source(paste(tdir,"/test_rom_o2.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_o3-0.r *****")
cat("\n")
source(paste(tdir,"/test_rom_o3-0.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_o3.r *****")
cat("\n")
source(paste(tdir,"/test_rom_o3.r",sep=""),local=TRUE)
cat("\n")

save.image(paste(tdir,"/Optical/.RData",sep=""))
