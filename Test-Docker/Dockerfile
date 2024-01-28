#
# Run desired global verification tests
#

# Use mpem/global image containing required
# R packages and global functions
FROM mpem/global

## make directory for test output
RUN mkdir -p /Test

## copy run file with desired tests
COPY ./tests.r Test/.

## source run file
CMD Rscript /Test/tests.r > \
            /Test/tests.out
