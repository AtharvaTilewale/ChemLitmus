"""
Logging configuration for ChemLitmus.

Sets up file and console logging with appropriate formatters and levels.
"""

import logging
import logging.handlers
from typing import Optional

from chemlitmus.config import settings


def setup_logging(
    name: Optional[str] = None,
    level: Optional[str] = None,
    log_file: Optional[str] = None,
) -> logging.Logger:
    """
    Configure logging for ChemLitmus.

    Args:
        name: Logger name (defaults to module name)
        level: Logging level (defaults to settings.log_level)
        log_file: Log file path (defaults to settings.log_file)

    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name or __name__)
    level = level or settings.log_level
    log_file = log_file or str(settings.log_file)

    # Avoid duplicate handlers
    if logger.handlers:
        return logger

    logger.setLevel(getattr(logging, level.upper()))

    # Formatter with timestamps and context
    formatter = logging.Formatter(
        fmt="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # File handler (rotating to avoid huge files)
    file_handler = logging.handlers.RotatingFileHandler(
        log_file,
        maxBytes=10 * 1024 * 1024,  # 10 MB
        backupCount=5,
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Console handler (less verbose)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.WARNING)
    console_formatter = logging.Formatter("%(levelname)s: %(message)s")
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    return logger


# Root logger for the package
logger = setup_logging("chemlitmus")
