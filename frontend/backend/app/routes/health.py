"""
Health check route.
"""

from __future__ import annotations

from fastapi import APIRouter

from app.config import settings

router = APIRouter()


@router.get("/health", tags=["System"])
async def health_check():
    """Verify that the API server is operational."""
    return {
        "status": "operational",
        "app_name": settings.APP_NAME,
        "version": settings.APP_VERSION,
    }
