"""MitoOmics-GPU package."""

__version__ = "0.2.0"

try:
    from .cli import main as cli_main  # noqa: F401
except Exception:
    def cli_main(*args, **kwargs):
        raise SystemExit("Use: python -m mitoomics_gpu ...")

try:
    from .mhi import combine_components  # noqa: F401
except Exception:
    pass

try:
    from .stats import mhi_group_test, mhi_correlation  # noqa: F401
except Exception:
    pass
