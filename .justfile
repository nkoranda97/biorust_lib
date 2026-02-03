dev:
    uv run maturin develop

check:
    cargo fmt --all
    uvx ruff format .
    uv run cargo clippy --all-targets -- -D warnings
    uv run cargo test --all
    uv run pytest -q