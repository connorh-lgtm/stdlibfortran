FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    gfortran g++ cmake ninja-build make python3 python3-pip git \
    && rm -rf /var/lib/apt/lists/*
ENV FC=gfortran FFLAGS="-O2 -fPIC" CFLAGS="-O2 -fPIC"
RUN python3 -m pip install --no-cache-dir pytest pytest-xdist fypp numpy joblib
WORKDIR /workspace
