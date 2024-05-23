## @brief  Dockerfile to containerise the IFM
##
FROM --platform=linux/amd64 ubuntu:22.04

RUN apt-get update
RUN apt-get install -y build-essential python3 python3-pip libnetcdf-dev libnetcdf-cxx-legacy-dev vim
RUN apt-get -y install netcdf-bin

RUN pip install numpy
RUN pip install netcdf4
RUN pip install utm
 
COPY . .

RUN make clean
RUN make -j

CMD ["ls", "-al"]
