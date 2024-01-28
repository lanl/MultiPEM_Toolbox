#
# Run desired Seismic sub-analysis
#

# Use iydt/seismic image containing required
# R packages, global functions, IYDT-specific
# functions and data, and seismic-specific
# functions
FROM iydt/seismic

## source run file
CMD Rscript /I-SUGAR-hob-pi-0/runMPEM.r > /I-SUGAR-hob-pi-0/runMPEM.out; \
Rscript /I-SUGAR-hob-pi-0/runMPEM_0.r > /I-SUGAR-hob-pi-0/runMPEM_0.out
