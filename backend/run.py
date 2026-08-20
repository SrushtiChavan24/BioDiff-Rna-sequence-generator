#!/usr/bin/env python3
"""
RNA Genomics Platform — Server Entry Point.

Usage:
    python run.py

Starts the FastAPI application via Uvicorn with settings drawn from
environment variables or app.config defaults.
"""

from __future__ import annotations

import uvicorn

from app.config import settings

if __name__ == "__main__":
    uvicorn.run(
        "app.main:app",
        host=settings.HOST,
        port=settings.PORT,
        reload=settings.DEBUG,
        log_level="info" if not settings.DEBUG else "debug",
    )
