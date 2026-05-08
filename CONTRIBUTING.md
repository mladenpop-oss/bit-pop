# Contributing to Bit-Pop

Thank you for your interest in contributing to Bit-Pop! This document covers how to get started, coding standards, and the contribution workflow.

## Getting Started

### Prerequisites

- Rust toolchain (2021 edition) — install via [rustup](https://rustup.rs)
- Python 3.x with Biopython — only for read simulation scripts
- Git

### Setup

```bash
git clone https://github.com/mladenpop-oss/bit-pop.git
cd bit-pop
cargo build --release
cargo test
```

All tests must pass before submitting a PR.

## Development Workflow

1. Create a feature branch: `git checkout -b feature/my-feature`
2. Make changes, write tests, run `cargo test` and `cargo clippy`
3. Commit with conventional commit messages: `feat:`, `fix:`, `docs:`, `refactor:`, `test:`, `chore:`
4. Push and open a Pull Request

## Coding Standards

- **No clippy warnings**: Run `cargo clippy -- -D warnings` before submitting
- **No build warnings**: Run `cargo build --release` and fix all warnings
- **Tests**: Unit tests for new functionality, integration tests for pipeline changes
- **2-bit DNA encoding**: Always use the 2-bit encoding throughout (A=0, C=1, G=2, T=3)
- **No unsafe**: Avoid `unsafe` blocks unless absolutely necessary and well-documented
- **Documentation**: Add doc comments to all public API items

## Adding New Features

1. Open an issue first to discuss the feature
2. Write tests first (TDD approach)
3. Implement the feature
4. Add benchmarks if the feature affects performance
5. Update README.md if the feature changes user-facing behavior

## Benchmarking

To run the benchmark suite:

```bash
cargo bench
```

When proposing performance improvements, include benchmark results comparing before/after.

## Testing

```bash
# All tests
cargo test

# Integration tests only
cargo test --test integration_tests

# Specific test
cargo test test_two_bit_align
```

Target: 80%+ test coverage.

## Code Review

PRs will be reviewed by project maintainers. Expect feedback on:

- Correctness of bioinformatics algorithms
- Performance implications
- Test coverage
- Code style and clippy compliance

## Questions?

Open a GitHub issue with the `question` label or start a GitHub Discussion.
