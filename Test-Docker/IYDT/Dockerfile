#
# Copy code for desired verification tests
#

# Use mpem/iydt image containing required
# R packages, global functions, and
# IYDT-specific functions and data
FROM mpem/iydt

## make directory for test output
RUN mkdir -p /Test

## copy run files for each test
COPY ./Applications/Test/IYDT/*.r Test/.
