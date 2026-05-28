#!/bin/bash
# Bit-Pop Bioconda build script
# Extracts pre-built binary from GitHub Release and installs to conda prefix

set -ex

# Extract the binary (already downloaded from source/url)
tar xzf *.tar.gz

# Install to conda bin directory
install -d $PREFIX/bin
install -m 755 bit-pop $PREFIX/bin/bit-pop
