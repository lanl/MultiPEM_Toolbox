#
# Create iydt/4-phen image.
#

# Use mpem/iydt image containing required
# R packages, global functions, and
# IYDT-specific functions and data
FROM mpem/iydt

# make directories for each sub-analysis
RUN mkdir -p /I-EIV
RUN mkdir -p /I-EIV-SUGAR-hob
RUN mkdir -p /I-EIV-SUGAR-hob-pi
RUN mkdir -p /4phen_eiv_cal
RUN mkdir -p /I-EIV-SUGAR-hob-0
RUN mkdir -p /I-EIV-SUGAR-hob-pi-0
# starting values: comment out if not
# available
RUN mkdir -p /Opt

# copy run files into each sub-analysis
COPY ./I-EIV/runMPEM.r I-EIV/.
COPY ./I-EIV-SUGAR-hob/runMPEM.r I-EIV-SUGAR-hob/.
COPY ./I-EIV-SUGAR-hob-pi/runMPEM.r I-EIV-SUGAR-hob-pi/.
#COPY ./I-EIV-SUGAR-hob-0/runMPEM.r I-EIV-SUGAR-hob-0/.
#COPY ./I-EIV-SUGAR-hob-pi-0/runMPEM.r I-EIV-SUGAR-hob-pi-0/.
COPY ./I-EIV-SUGAR-hob-0/runMPEM.r 4phen_eiv_cal/.
COPY ./I-EIV-SUGAR-hob-0/runMPEM_0.r I-EIV-SUGAR-hob-0/.
COPY ./I-EIV-SUGAR-hob-pi-0/runMPEM_0.r I-EIV-SUGAR-hob-pi-0/.
# starting values: comment out if not
# available
COPY ./Opt/opt_1.RData Opt/.
COPY ./Opt/opt_2.RData Opt/.
COPY ./Opt/opt_3_eiv.RData Opt/.
COPY ./Opt/opt_4_eiv.RData Opt/.
COPY ./Opt/opt_1_0.RData Opt/.
COPY ./Opt/opt_2_0.RData Opt/.
COPY ./Opt/opt_3_eiv_0.RData Opt/.
COPY ./Opt/opt_4_eiv_0.RData Opt/.
