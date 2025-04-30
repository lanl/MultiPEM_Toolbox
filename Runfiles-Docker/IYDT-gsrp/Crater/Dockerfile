#
# Create iydt/crater image.
#

# Use mpem/iydt image containing required
# R packages, global functions, and
# IYDT-specific functions and data
FROM mpem/iydt

# make directories for each sub-analysis
RUN mkdir -p /I-EIV-SUGAR
RUN mkdir -p /crater_eiv_cal
RUN mkdir -p /I-EIV-SUGAR-0
RUN mkdir -p /I-SUGAR
RUN mkdir -p /crater_cal
RUN mkdir -p /I-SUGAR-0
# starting values: comment out if not
# available
#RUN mkdir -p /Opt

# copy run files into each sub-analysis
COPY ./I-EIV-SUGAR/runMPEM.r I-EIV-SUGAR/.
#COPY ./I-EIV-SUGAR-0/runMPEM.r I-EIV-SUGAR-0/.
COPY ./I-EIV-SUGAR-0/runMPEM.r crater_eiv_cal/.
COPY ./I-EIV-SUGAR-0/runMPEM_0.r I-EIV-SUGAR-0/.
COPY ./I-SUGAR/runMPEM.r I-SUGAR/.
#COPY ./I-SUGAR-0/runMPEM.r I-SUGAR-0/.
COPY ./I-SUGAR-0/runMPEM.r crater_cal/.
COPY ./I-SUGAR-0/runMPEM_0.r I-SUGAR-0/.
# starting values: comment out if not
# available
#COPY ./Opt/opt.RData Opt/.
