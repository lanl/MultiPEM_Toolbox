#
# Run desired Crater sub-analysis
#

# Use iydt/crater image containing required
# R packages, global functions, IYDT-specific
# functions and data, and crater-specific
# functions
FROM iydt/crater

## source run files
CMD Rscript /I-EIV-SUGAR-0/runMPEM.r > /I-EIV-SUGAR-0/runMPEM.out; \
Rscript /I-EIV-SUGAR-0/runMPEM_0.r > /I-EIV-SUGAR-0/runMPEM_0.out
