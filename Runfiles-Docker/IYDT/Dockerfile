#
# Create mpem/iydt image
#

# Use mpem/global image containing required
# R packages and global functions
FROM mpem/global

## create directories
RUN mkdir -p /Code/IYDT
RUN mkdir -p /Data

# copy IYDT files
COPY ./Applications/Code/IYDT/*.r Code/IYDT/.
COPY ./Applications/Code/IYDT/Phenomenology/*.r Code/IYDT/.
COPY ./Applications/Data/IYDT/*.csv Data/.
