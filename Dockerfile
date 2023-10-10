FROM ubuntu:20.04
LABEL maintainer="Schuyler <schuyler.smith@nutrien.com>"


ENV DEBIAN_FRONTEND noninteractive
RUN apt update \
    && apt install --yes --no-install-recommends \
        ca-certificates \
        wget \
        git \
        zip unzip\
        libbz2-dev \
        libpq-dev \
        g++ \
        make \
        cmake \
        libboost-all-dev \
    && apt-get upgrade --yes \
    && apt-get clean

WORKDIR /opt
ENV PATH="${PATH}:/opt/"


CMD ["/bin/bash"]
