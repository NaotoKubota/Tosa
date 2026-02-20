###############################################################
<<<<<<< HEAD
# Dockerfile to build container images for Tosa
=======
# Dockerfile to build container images for Tosa (v0.3.0)
>>>>>>> 33a432e75c1c3be5acd1d209fd1f1af23b49563d
# Based on rust slim-bullseye
################################################################

FROM rust:slim-bullseye

# File Author / Maintainer
LABEL maintainer="Naoto Kubota <naotok@ucr.edu>"

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies first
RUN apt-get -qq update && \
	apt-get -qq -y install \
    build-essential gcc g++ make cmake libtool texinfo dpkg-dev pkg-config \
    libgsl-dev wget locate less vim zlib1g-dev bzip2 lzma curl perl \
    libboost-dev libcurl4-openssl-dev libboost-all-dev libbz2-dev liblzma-dev \
    libpcre3 libpcre3-dev libclang-dev

<<<<<<< HEAD
# Copy source code and build
COPY . /opt_tosa/Tosa
RUN cd /opt_tosa/Tosa && \
=======
# Install Tosa
RUN mkdir -p /opt_tosa && \
    cd /opt_tosa && \
    git clone https://github.com/NaotoKubota/Tosa.git -b v0.3.0 && \
    cd Tosa && \
>>>>>>> 33a432e75c1c3be5acd1d209fd1f1af23b49563d
    cargo build --release && \
    cp target/release/tosa /usr/local/bin/ && \
    rm -rf /opt_tosa

# Clean up
RUN apt-get -qq -y remove build-essential gcc g++ make cmake libtool texinfo dpkg-dev pkg-config && \
    apt-get -qq -y autoremove && \
    apt-get -qq -y clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set working directory
WORKDIR /home

# bash
CMD ["bash"]
