# .devcontainer/Dockerfile

# Use the official Ubuntu as a parent image
FROM mcr.microsoft.com/vscode/devcontainers/base:ubuntu-20.04

# Prevent interactive prompts during package installation
ARG DEBIAN_FRONTEND=noninteractive

# Update package lists and install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    g++-10 \
    clang \
    pkg-config \    
    && rm -rf /var/lib/apt/lists/*

FROM hdfgroup/hdf5lib:latest

# Set g++-10 as the default compiler if available
RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100