#
# Run desired Optical sub-analysis
#

# Use iydt/optic image containing required
# R packages, global functions, IYDT-specific
# functions and data, and optical-specific
# functions
FROM iydt/optic

## source run files
CMD Rscript /I-EIV-SUGAR-hob-pi-0/runMPEM.r > /I-EIV-SUGAR-hob-pi-0/runMPEM.out; \
Rscript /I-EIV-SUGAR-hob-pi-0/runMPEM_0.r > /I-EIV-SUGAR-hob-pi-0/runMPEM_0.out
