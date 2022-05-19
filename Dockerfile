FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=US/Central

COPY . /bigstick
COPY clean_src /bigstick/src

WORKDIR /bigstick/src

# Install dependencies
RUN apt-get update
RUN apt-get install -y make cmake gfortran nano python3 python3-pip
RUN apt-get clean -q
RUN pip3 install plotly streamlit


# build and compile bigstick
RUN make gfortran
WORKDIR /bigstick