"""Configuration loading.

The single entry point for analysis parameters. Stages take a :class:`Config`
and never read ``configs/analysis.yaml`` themselves, so a test can hand a stage
a modified config without touching the file on disk.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


def project_root() -> Path:
    """Repository root, resolved from this file's location."""
    return Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Config:
    """Parsed ``analysis.yaml`` with path helpers bound to the repository root."""

    data: dict[str, Any]
    root: Path
    source: Path | None = None

    # -- dict-ish access -------------------------------------------------
    def __getitem__(self, key: str) -> Any:
        return self.data[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self.data.get(key, default)

    # -- path helpers ----------------------------------------------------
    def path(self, kind: str, *parts: str) -> Path:
        """Absolute path inside one of the directories declared under ``paths``.

        ``cfg.path("results", "markers.csv")`` -> ``<root>/results/markers.csv``.
        The parent directory is created so callers can write immediately.
        """
        try:
            base = self.data["paths"][kind]
        except KeyError as exc:  # pragma: no cover - configuration error
            raise KeyError(
                f"unknown path kind {kind!r}; config declares {sorted(self.data['paths'])}"
            ) from exc
        out = self.root / base
        out.mkdir(parents=True, exist_ok=True)
        return out.joinpath(*parts) if parts else out

    @property
    def seed(self) -> int:
        return int(self.data.get("seed", 0))


def load_config(path: str | os.PathLike[str] | None = None) -> Config:
    """Load ``configs/analysis.yaml`` (or an explicit path) into a :class:`Config`."""
    root = project_root()
    src = Path(path) if path is not None else root / "configs" / "analysis.yaml"
    if not src.is_absolute():
        src = (root / src).resolve()
    with open(src, "r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"{src} did not parse to a mapping")
    return Config(data=data, root=root, source=src)
