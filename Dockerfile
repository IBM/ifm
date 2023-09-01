## @brief  Dockerfile to containerise shared-memory version of IFM
##
## @author Maksims Abalenkovs
## @email  maksims.abalenkovs@stfc.ac.uk
## @date   Sep 1, 2023
## @version 0.1

FROM ubuntu:22.04
RUN apt-get update
RUN apt-get install -y build-essential git libnetcdf-dev libnetcdf-cxx-legacy-dev
ARG GIT_USERNAME
ARG GIT_ACCESS_TOKEN
RUN git clone https://${GIT_USERNAME}:${GIT_ACCESS_TOKEN}@gitlab.stfc.ac.uk/cimf/ifm.git /ifm/
WORKDIR /ifm
RUN make -j

## @todo Replace with actual call to IFM executable, e.g. CMD ["./app.py"]
CMD ["ls", "-al"]

## @eof Dockerfile
