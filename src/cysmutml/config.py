"""Configuration loading."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_config(path: str | Path = "configs/default.yaml") -> dict[str, Any]:
    with Path(path).open() as handle:
        return yaml.safe_load(handle)
