########################################################################
#                                                                      #
##' Catch *and* save both errors and warnings, and in the case of      #
##' a warning, also keep the computed result.                          #
##'                                                                    #
##' @title tryCatch both warnings (with value) and errors              #
##' @param expr an \R expression to evaluate                           #
##' @return a list with 'value' and 'warning', where                   #
##'   'value' may be an error caught.                                  #
##' @author Martin Maechler;                                           #
##' Copyright (C) 2010-2012  The R Core Team                           #
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

tryCatch.W.E <- function(expr)
{
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value=withCallingHandlers(tryCatch(expr,error=function(e) e),
                                 warning=w.handler), warning=W)
}
