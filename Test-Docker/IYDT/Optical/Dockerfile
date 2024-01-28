#
# Run desired Optical verification tests
#

# Use iydt/test image containing required
# functions for verification testing
FROM iydt/test

## make directory for test output
RUN mkdir -p /Test/Optical

## copy run file with desired tests
COPY ./tests.r Test/Optical/.

## source run file
CMD Rscript /Test/Optical/tests.r > \
            /Test/Optical/tests.out
