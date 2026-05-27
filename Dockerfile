# Bit-Pop Dockerfile
# Multi-stage build: Rust compilation → minimal runtime image

# ============================================
# Stage 1: Build
# ============================================
FROM rust:1.87-slim AS builder

WORKDIR /usr/src/bit-pop

# Install build dependencies (libsais C library, etc.)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        pkg-config \
        libz-dev && \
    rm -rf /var/lib/apt/lists/*

# Copy all source code
COPY . .

# Build release binary
RUN cargo build --release && \
    cp target/release/bit-pop /usr/local/bin/bit-pop

# ============================================
# Stage 2: Runtime
# ============================================
FROM debian:bookworm-slim AS runtime

LABEL org.opencontainers.image.title="Bit-Pop" \
      org.opencontainers.image.description="Ultra-fast multi-genome DNA read classification" \
      org.opencontainers.image.url="https://github.com/mladenpop-oss/bit-pop" \
      org.opencontainers.image.source="https://github.com/mladenpop-oss/bit-pop" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.vendor="Mladen Popovic"

# Install minimal runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        libssl3 && \
    rm -rf /var/lib/apt/lists/* && \
    useradd --create-home --shell /bin/bash bitpop

# Copy binary from builder
COPY --from=builder /usr/local/bin/bit-pop /usr/local/bin/bit-pop

# Set working directory and permissions
WORKDIR /home/bitpop
RUN chown -R bitpop:bitpop /home/bitpop

USER bitpop

# Default command
ENTRYPOINT ["bit-pop"]
CMD ["--help"]
