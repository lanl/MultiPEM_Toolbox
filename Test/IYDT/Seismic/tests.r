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
adir = "../Code"

# source forward models and jacobians
source(paste(adir,"/forward.r",sep=""),local=TRUE)
source(paste(adir,"/forward_0.r",sep=""),local=TRUE)
source(paste(adir,"/jacobian.r",sep=""),local=TRUE)
source(paste(adir,"/jacobian_0.r",sep=""),local=TRUE)

# directory to test files
tdir = "../Test"

print("***** test_rom_s1_d.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s1_d.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s1_v.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s1_v.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s1a_d.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s1a_d.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s1a_v.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s1a_v.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s2_d-0.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s2_d-0.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s2_d.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s2_d.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s2_v-0.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s2_v-0.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s2_v.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s2_v.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s2a_d-0.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s2a_d-0.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s2a_d.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s2a_d.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s2a_v-0.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s2a_v-0.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s2a_v.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s2a_v.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s2b_d-0.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s2b_d-0.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s2b_d.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s2b_d.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s2b_v-0.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s2b_v-0.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s2b_v-0.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s2b_v-0.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s3_d-0.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s3_d-0.r",sep=""),local=TRUE)
cat("\n")

print("***** test_rom_s3_d.r *****")
cat("\n")
source(paste(tdir,"/test_rom_s3_d.r",sep=""),local=TRUE)
