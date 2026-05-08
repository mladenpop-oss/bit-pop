#!/bin/bash

if [ "$(uname)" == "Darwin" ]; then
  # macOS - use the macOS binary
  tar xzf bit-pop-x86_64-macos.tar.gz 2>/dev/null || \
  tar xzf bit-pop-aarch64-macos.tar.gz
  install -m 755 bit-pop $PREFIX/bin/bit-pop
else
  # Linux
  tar xzf bit-pop-x86_64-linux.tar.gz
  install -m 755 bit-pop $PREFIX/bin/bit-pop
fi
