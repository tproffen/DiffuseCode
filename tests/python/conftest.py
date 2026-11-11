"""Register the manifest tags as pytest markers (tests/manifest.toml)."""

TAGS = {
    "fast": "runs in under 10 s",
    "medium": "runs in under 30 s",
    "slow": "runs for 30 s or more",
    "mc": "stochastic (Monte Carlo / random); results vary between runs",
    "mpi": "needs MPI",
    "nocheck": "has no self-check yet",
    "exclusive": "sets its own OpenMP thread count; do not run in parallel",
}


def pytest_configure(config):
    for tag, text in TAGS.items():
        config.addinivalue_line("markers", f"{tag}: {text}")
