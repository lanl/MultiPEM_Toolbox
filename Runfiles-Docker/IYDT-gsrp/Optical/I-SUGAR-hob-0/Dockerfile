#
# Run desired Optical sub-analysis
#

# Use iydt/optic image containing required
# R packages, global functions, IYDT-specific
# functions and data, and optical-specific
# functions
FROM iydt/optic

## source run files
CMD Rscript /I-SUGAR-hob-0/runMPEM.r > /I-SUGAR-hob-0/runMPEM.out; \
Rscript /I-SUGAR-hob-0/runMPEM_0.r > /I-SUGAR-hob-0/runMPEM_0.out
