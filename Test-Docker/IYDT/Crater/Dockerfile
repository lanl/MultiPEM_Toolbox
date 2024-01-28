#
# Run desired Crater verification tests
#

# Use iydt/test image containing required
# functions for verification testing
FROM iydt/test

## make directory for test output
RUN mkdir -p /Test/Crater

## copy run file with desired tests
COPY ./tests.r Test/Crater/.

## source run file
CMD Rscript /Test/Crater/tests.r > \
            /Test/Crater/tests.out
