FROM ubuntu:18.04
# maybe try alpine. ubuntu is a bit heavy

COPY src_bigstick /src_bigstick

WORKDIR /src_bigstick

# Install dependencies
RUN apt-get update
RUN apt-get install -y make gfortran
RUN apt-get clean -q


# build and compile bigstick
RUN make gfortran
# RUN ls | grep bigstick