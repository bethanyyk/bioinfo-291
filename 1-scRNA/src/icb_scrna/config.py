"""Configuration loading.

Uses the reproducible-skills Config base class (config_base.py) for parameter
management. Stages take a :class:`Config` and never read ``configs/analysis.yaml``
themselves, so a test can hand a stage a modified config without touching the file
on disk.

See ../config_base.py (from reproducible-skills) for the base implementation.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

# Import base Config class from reproducible-skills (copied into src/)
from config_base import Config as BaseConfig, load_config as base_load


def project_root() -> Path:
    """Repository root, resolved from this file's location."""
    return Path(__file__).resolve().parents[2]


class Config(BaseConfig):
    """scRNA-specific configuration wrapping reproducible-skills base.
    
    Extends BaseConfig with domain-specific convenience methods while keeping
    parameter loading, path resolution, and dict-like access from the base class.
    """

    @property
    def seed(self) -> int:
        """Random seed for reproducibility (default: 0)."""
        return int(self.data.get("seed", 0))


def load_config(path: str | os.PathLike[str] | None = None) -> Config:
    """Load ``configs/analysis.yaml`` (or an explicit path) into a :class:`Config`.
    
    Uses reproducible-skills base loader, wrapping result in scRNA Config class.
    """
    root = project_root()
    src = Path(path) if path is not None else root / "configs" / "analysis.yaml"
    if not src.is_absolute():
        src = (root / src).resolve()
    
    # Load using base class
    base_cfg = base_load(src)
    
    # Wrap in scRNA Config to add domain-specific properties
    cfg = Config(data=base_cfg.data, root=base_cfg.root, source=src)
    return cfg
