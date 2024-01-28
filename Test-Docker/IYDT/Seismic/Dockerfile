#
# Run desired Seismic verification tests
#

# Use iydt/test image containing required
# functions for verification testing
FROM iydt/test

## make directory for test output
RUN mkdir -p /Test/Seismic

## copy run file with desired tests
COPY ./tests.r Test/Seismic/.

## source run file
CMD Rscript /Test/Seismic/tests.r > \
            /Test/Seismic/tests.out
