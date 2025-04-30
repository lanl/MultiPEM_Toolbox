#
# Create mpem/iydt image
#

# Use mpem/global image containing required
# R packages and global functions
FROM mpem/global

## create directories
RUN mkdir -p /Code/IYDT-gsrp
RUN mkdir -p /Data

# copy IYDT files
COPY ./Applications/Code/IYDT-gsrp/*.r Code/IYDT-gsrp/.
COPY ./Applications/Code/IYDT-gsrp/Phenomenology/*.r Code/IYDT-gsrp/.
COPY ./Applications/Data/IYDT-gsrp/*.csv Data/.
