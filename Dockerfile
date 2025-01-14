FROM ubuntu:latest

# Install required packages
RUN apt-get update && apt-get install -y \
    openjdk-11-jdk \
    curl \
    wget \
    git \
    procps

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash
RUN mv nextflow /usr/local/bin/

# Set working directory
WORKDIR /workspace