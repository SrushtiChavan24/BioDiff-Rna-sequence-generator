from __future__ import annotations

from datetime import datetime
from typing import Any, Generic, TypeVar

from pydantic import BaseModel, Field

T = TypeVar("T")


class BaseResponse(BaseModel, Generic[T]):
    """Standard API response envelope for all endpoints.

    Provides a consistent response structure across the entire API:
    - success: whether the request was processed without error
    - message: human-readable summary
    - data: the actual response payload (generic type)
    - timestamp: ISO 8601 time of response generation
    - version: API version identifier
    """

    success: bool = Field(
        default=True,
        description="Indicates whether the request completed successfully.",
    )
    message: str = Field(
        default="Request completed successfully.",
        description="Human-readable status message.",
    )
    data: T | None = Field(
        default=None,
        description="Response payload, typed per-endpoint.",
    )
    timestamp: str = Field(
        default_factory=lambda: datetime.utcnow().isoformat() + "Z",
        description="ISO 8601 UTC timestamp of response generation.",
    )
    version: str = Field(
        default="1.0.0",
        description="API version identifier.",
    )


class ErrorResponse(BaseModel):
    """Standard error response for API failures."""

    success: bool = Field(
        default=False,
        description="Always false for error responses.",
    )
    detail: str = Field(
        ...,
        description="Human-readable error description.",
    )
    status_code: int = Field(
        ...,
        description="HTTP status code of the error.",
    )
    timestamp: str = Field(
        default_factory=lambda: datetime.utcnow().isoformat() + "Z",
        description="ISO 8601 UTC timestamp.",
    )
