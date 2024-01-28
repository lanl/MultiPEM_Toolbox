#
# Create iydt/seismic image.
#

# Use mpem/iydt image containing required
# R packages, global functions, and
# IYDT-specific functions and data
FROM mpem/iydt

# make directories for each sub-analysis
RUN mkdir -p /I
RUN mkdir -p /I-SUGAR-hob
RUN mkdir -p /I-SUGAR-hob-pi
RUN mkdir -p /seismic_cal
RUN mkdir -p /I-SUGAR-hob-0
RUN mkdir -p /I-SUGAR-hob-pi-0
# starting values: comment out if not
# available
#RUN mkdir -p /Opt

# copy run files into each sub-analysis
COPY ./I/runMPEM.r I/.
COPY ./I-SUGAR-hob/runMPEM.r I-SUGAR-hob/.
COPY ./I-SUGAR-hob-pi/runMPEM.r I-SUGAR-hob-pi/.
#COPY ./I-SUGAR-hob-0/runMPEM.r I-SUGAR-hob-0/.
#COPY ./I-SUGAR-hob-pi-0/runMPEM.r I-SUGAR-hob-pi-0/.
COPY ./I-SUGAR-hob-0/runMPEM.r seismic_cal/.
COPY ./I-SUGAR-hob-0/runMPEM_0.r I-SUGAR-hob-0/.
COPY ./I-SUGAR-hob-pi-0/runMPEM_0.r I-SUGAR-hob-pi-0/.
# starting values: comment out if not
# available
#COPY ./Opt/opt.RData Opt/.
