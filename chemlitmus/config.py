"""
Configuration management for ChemLitmus.

Handles environment variables, cache directories, API endpoints, and defaults.
"""

import os
from pathlib import Path

import platformdirs
from pydantic import ConfigDict, Field
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """
    Application settings with environment variable support.

    Environment variables override defaults. Example:
        export CHEMLITMUS_CACHE_DIR=/custom/cache
        export CHEMLITMUS_LOG_LEVEL=DEBUG
    """

    # Directories
    cache_dir: Path = Field(
        default_factory=lambda: Path(platformdirs.user_cache_dir("chemlitmus")),
        description="Cache directory for database and downloads",
    )
    data_dir: Path = Field(
        default_factory=lambda: Path(platformdirs.user_data_dir("chemlitmus")),
        description="Data directory for user files",
    )
    log_dir: Path = Field(
        default_factory=lambda: Path(platformdirs.user_log_dir("chemlitmus")),
        description="Log directory",
    )

    # API settings
    pubchem_base_url: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    pubchem_timeout: int = Field(default=10, description="HTTP timeout in seconds")
    pubchem_retries: int = Field(default=3, description="Number of retry attempts")

    # Performance
    batch_size: int = Field(default=50, description="Batch size for bulk operations")
    max_workers: int = Field(
        default_factory=os.cpu_count,
        description="Maximum number of concurrent threads/workers",
    )
    rate_limit_delay: float = Field(
        default=0.5,
        description="Delay between requests to PubChem (seconds)",
    )

    # Database
    db_name: str = "chemlitmus.db"
    enable_cache: bool = True

    # Logging
    log_level: str = Field(default="INFO", description="Logging level")
    log_format: str = (
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Feature flags
    enable_async: bool = True
    enable_rdkit_validation: bool = True

    model_config = ConfigDict(
        env_prefix="CHEMLITMUS_",
        case_sensitive=False,
        env_file=".env",
    )

    def __init__(self, **kwargs):
        """Initialize settings and create directories."""
        super().__init__(**kwargs)
        self._create_directories()

    def _create_directories(self) -> None:
        """Create required directories if they don't exist."""
        for directory in [self.cache_dir, self.data_dir, self.log_dir]:
            directory.mkdir(parents=True, exist_ok=True)

    @property
    def db_path(self) -> Path:
        """Path to SQLite database."""
        return self.cache_dir / self.db_name

    @property
    def log_file(self) -> Path:
        """Path to main log file."""
        return self.log_dir / "chemlitmus.log"


# Global settings instance
settings = Settings()
