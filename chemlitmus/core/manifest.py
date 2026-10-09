"""Run manifests: what produced a result, precisely enough to reproduce or to detect drift."""

from __future__ import annotations

import hashlib
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

from pydantic import BaseModel, Field

MANIFEST_SCHEMA_VERSION = "1"


def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


class FileRef(BaseModel):
    path: str = Field(description="File name only (no private directory structure).")
    sha256: str
    size_bytes: int


class RunManifest(BaseModel):
    schema_version: str = MANIFEST_SCHEMA_VERSION
    tool: str = "chemlitmus"
    tool_version: str
    command: str
    created_at: str
    python_version: str
    rdkit_version: Optional[str] = None
    platform: str
    inputs: List[FileRef] = Field(default_factory=list)
    outputs: List[FileRef] = Field(default_factory=list)
    policy_hash: Optional[str] = None
    policy_name: Optional[str] = None
    rule_catalogue: Optional[Dict[str, str]] = Field(None, description="identity/hash of the alert catalogue used")
    reference_library: Optional[Dict[str, str]] = None
    seeds: Dict[str, int] = Field(default_factory=dict)
    settings: Dict[str, object] = Field(default_factory=dict, description="Runtime settings relevant to the result (no secrets, no private paths).")
    output_schema_version: Optional[str] = None
    notes: List[str] = Field(default_factory=list)

    def add_input(self, path: str | Path) -> None:
        p = Path(path)
        self.inputs.append(FileRef(path=p.name, sha256=sha256_file(p), size_bytes=p.stat().st_size))

    def add_output(self, path: str | Path) -> None:
        p = Path(path)
        if p.exists():
            self.outputs.append(FileRef(path=p.name, sha256=sha256_file(p), size_bytes=p.stat().st_size))

    def save(self, path: str | Path) -> None:
        Path(path).write_text(self.model_dump_json(indent=2) + "\n")


def new_manifest(command: str, **settings) -> RunManifest:
    from chemlitmus import __version__
    try:
        import rdkit
        rd = rdkit.__version__
    except Exception:  # pragma: no cover
        rd = None
    return RunManifest(
        tool_version=__version__, command=command, created_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        python_version=sys.version.split()[0], rdkit_version=rd, platform=f"{platform.system()} {platform.machine()}", settings=settings,
    )


__all__ = ["RunManifest", "FileRef", "new_manifest", "sha256_file", "sha256_text", "MANIFEST_SCHEMA_VERSION"]
