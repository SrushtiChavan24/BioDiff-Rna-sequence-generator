from __future__ import annotations

import os
from pathlib import Path

from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.orm import DeclarativeBase

from app.config import settings

# ---------------------------------------------------------------------------
# Async SQLite Engine & Session Factory
# ---------------------------------------------------------------------------
# Uses aiosqlite as the async driver for SQLite.  The database file lives
# under the project-level `outputs/` directory so it persists across restarts
# and can be included in exports.

DB_DIR: Path = settings.OUTPUTS_DIR
DB_DIR.mkdir(parents=True, exist_ok=True)

DATABASE_URL: str = os.getenv(
    "DATABASE_URL",
    f"sqlite+aiosqlite:///{DB_DIR / 'rna_genomics.db'}",
)

_engine = create_async_engine(
    DATABASE_URL,
    echo=settings.DEBUG,
    connect_args={"check_same_thread": False},
)

_async_session_factory = async_sessionmaker(
    bind=_engine,
    class_=AsyncSession,
    expire_on_commit=False,
)


# ---------------------------------------------------------------------------
# Declarative Base
# ---------------------------------------------------------------------------

class Base(DeclarativeBase):
    """SQLAlchemy declarative base for all ORM models."""


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

async def get_async_session() -> AsyncSession:  # type: ignore[misc]
    """Yield an async SQLAlchemy session — intended for FastAPI dependency injection."""
    async with _async_session_factory() as session:
        try:
            yield session
            await session.commit()
        except Exception:
            await session.rollback()
            raise
        finally:
            await session.close()


async def init_db() -> None:
    """Create all tables defined in the ORM metadata.

    Idempotent — safe to call on every startup.
    """
    async with _engine.begin() as conn:
        from app.database.models import AnalysisModel  # noqa: F401 – ensure models are loaded
        await conn.run_sync(Base.metadata.create_all)


async def close_db() -> None:
    """Dispose of the engine — call during application shutdown."""
    await _engine.dispose()
