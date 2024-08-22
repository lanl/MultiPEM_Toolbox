########################################################################
#                                                                      #
# Compute Jacobian from Cholesky factor elements to Covariance         #
# elements                                                             #
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

require(Matrix)

set.seed(110)

# Transform function
t_chol <- function(l, R)
{
  L = matrix(0,R,R)
  L[upper.tri(L,diag=TRUE)] = l
  diag(L) = exp(diag(L))

  Sigma = t(L) %*% L
  varp = diag(Sigma)
  sdp = sqrt(varp)
  
  tc = numeric(R*(R+1)/2)
  qq = 0
  for( r2 in 1:R ){
    for( r1 in 1:r2 ){
      qq = qq+1
      if( r1 == r2 ){
        tc[qq] = sum(L[1:r2,r2]^2)
      } else {
        tc[qq] = sum(L[1:r1,r1]*L[1:r1,r2])/sdp[r1]/
                 sdp[r2]
      }
    }
  }
  tc
}

vR = 1:10
for( R in vR ){
  npar = R*(R+1)/2

  l = rnorm(npar,mean=0,sd=3)

  L = Diagonal(R,exp(l[1:R]))
  if( R > 1 ){
    L[upper.tri(L)] = l[-(1:R)]
  }

  Sigma = t(L) %*% L
  varp = diag(Sigma)
  sdp = sqrt(varp)
  iS = Diagonal(R,1/sdp)
  C = iS %*% Sigma %*% iS

  k_ij = function(i,j)
  {
    i+choose(j,2)
  }

  Jac_c = Matrix(0,npar,npar,sparse=FALSE,doDiag=FALSE)

  qq = 0
  for( r2 in 1:R ){
    for( r1 in 1:r2 ){
      qq = qq+1
      if( r1 == r2 ){
        Jac_c[qq,k_ij(1,r2):k_ij(r2,r2)] = 2*L[1:r2,r2]
        Jac_c[qq,k_ij(r2,r2)] = Jac_c[qq,k_ij(r2,r2)]*L[r2,r2]
      } else {
        Jac_c[qq,k_ij(1,r1):k_ij(r1,r1)] =
        L[1:r1,r2]/sdp[r1]/sdp[r2]-C[r1,r2]*L[1:r1,r1]/varp[r1]
        Jac_c[qq,k_ij(r1,r1)] = Jac_c[qq,k_ij(r1,r1)]*L[r1,r1]
        Jac_c[qq,k_ij(1,r2):k_ij(r2,r2)] =
        c(L[1:r1,r1]/sdp[r1]/sdp[r2]-C[r1,r2]*L[1:r1,r2]/varp[r2],
          -C[r1,r2]*L[(r1+1):r2,r2]/varp[r2])
        Jac_c[qq,k_ij(r2,r2)] = Jac_c[qq,k_ij(r2,r2)]*L[r2,r2]
      }
    }
  }

  # Numerically check Jacobian

  tL = matrix(0,R,R)
  diag(tL) = l[1:R]
  if( R > 1 ){
    tL[upper.tri(tL)] = l[-(1:R)]
  }
  tl = tL[upper.tri(tL,diag=TRUE)]

  require(numDeriv)
  Jac_c_num = jacobian(t_chol, tl, method.args=list(r=4),
                       R=R) # numerical jacobian

  diff = Jac_c - Jac_c_num
  print(c(min(diff),max(diff)))
}
