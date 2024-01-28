#
# Run desired Acoustic verification tests
#

# Use iydt/test image containing required
# functions for verification testing
FROM iydt/test

## make directory for test output
RUN mkdir -p /Test/Acoustic

## copy run file with desired tests
COPY ./tests.r Test/Acoustic/.

## source run file
CMD Rscript /Test/Acoustic/tests.r > \
            /Test/Acoustic/tests.out
