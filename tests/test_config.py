"""
Tests for configuration module.

Verifies that settings are loaded correctly and directories are created.
"""

from pathlib import Path


from chemlitmus.config import Settings, settings


class TestSettings:
    """Test suite for Settings class."""

    def test_settings_instantiation(self, temp_cache_dir: Path) -> None:
        """Test that settings can be instantiated with custom paths."""
        test_settings = Settings(cache_dir=temp_cache_dir)
        assert test_settings.cache_dir == temp_cache_dir
        assert test_settings.cache_dir.exists()

    def test_environment_variable_override(self, temp_cache_dir: Path, monkeypatch) -> None:
        """Test that environment variables override defaults."""
        monkeypatch.setenv("CHEMLITMUS_LOG_LEVEL", "DEBUG")
        test_settings = Settings(cache_dir=temp_cache_dir)
        assert test_settings.log_level == "DEBUG"

    def test_db_path_property(self, test_settings: Settings) -> None:
        """Test that db_path property returns correct path."""
        expected = test_settings.cache_dir / test_settings.db_name
        assert test_settings.db_path == expected

    def test_log_file_property(self, test_settings: Settings) -> None:
        """Test that log_file property returns correct path."""
        expected = test_settings.log_dir / "chemlitmus.log"
        assert test_settings.log_file == expected

    def test_directories_created(self, test_settings: Settings) -> None:
        """Test that all required directories are created on init."""
        assert test_settings.cache_dir.exists()
        assert test_settings.data_dir.exists()
        assert test_settings.log_dir.exists()

    def test_pubchem_settings(self, test_settings: Settings) -> None:
        """Test PubChem configuration."""
        assert "pubchem.ncbi.nlm.nih.gov" in test_settings.pubchem_base_url
        assert test_settings.pubchem_timeout > 0
        assert test_settings.pubchem_retries > 0


class TestGlobalSettings:
    """Test suite for global settings instance."""

    def test_global_settings_exists(self) -> None:
        """Test that global settings instance is available."""
        assert settings is not None
        assert isinstance(settings, Settings)

    def test_global_settings_db_path(self) -> None:
        """Test that global settings has valid database path."""
        assert settings.db_path is not None
        assert settings.db_path.parent.exists()
