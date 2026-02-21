"""Logging configuration for OptiType."""

import logging
import time

# Package-level logger
logger = logging.getLogger("optitype")

# Record start time at import for elapsed-time display
_start_time = time.monotonic()


def elapsed() -> str:
    """Return elapsed time since package load as H:MM:SS.ss."""
    delta = time.monotonic() - _start_time
    minutes, seconds = divmod(delta, 60)
    hours, minutes = divmod(int(minutes), 60)
    return f"{hours}:{minutes:02d}:{seconds:05.2f}"


def configure_logging(verbose: bool) -> None:
    """Configure the optitype logger for CLI or API use."""
    level = logging.DEBUG if verbose else logging.WARNING
    logger.setLevel(level)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(handler)
