"""
Application configuration module.

Centralizes all environment-driven settings via Pydantic BaseSettings.
"""

from __future__ import annotations

import os
from pathlib import Path


class Settings:
    """Application settings populated from environment variables with sane defaults."""

    APP_NAME: str = "RNA Genomics Intelligence Platform"
    APP_VERSION: str = "1.0.0"
    APP_DESCRIPTION: str = (
        "Production-grade AI-powered RNA genomics analysis platform for "
        "therapeutic candidate ranking, structural analysis, and molecular "
        "similarity computation."
    )
    DEBUG: bool = os.getenv("DEBUG", "false").lower() == "true"
    API_PREFIX: str = "/api"

    # CORS
    CORS_ORIGINS: list[str] = [
        origin.strip()
        for origin in os.getenv(
            "CORS_ORIGINS", "http://localhost:3000"
        ).split(",")
    ]

    # Server
    HOST: str = os.getenv("HOST", "0.0.0.0")
    PORT: int = int(os.getenv("PORT", "8000"))

    # Analysis constraints
    MIN_SEQUENCE_LENGTH: int = 10
    MAX_SEQUENCE_LENGTH: int = 10_000
    MAX_BATCH_SIZE: int = 50

    # Paths
    BASE_DIR: Path = Path(__file__).resolve().parent.parent
    OUTPUTS_DIR: Path = BASE_DIR.parent / "outputs"

    def __init__(self) -> None:
        self.OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)


settings = Settings()
