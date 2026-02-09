"""Clustal Omega MSA integration."""

from __future__ import annotations

import importlib.resources
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def find_clustalo() -> Path:
    """Locate the clustalo binary.

    Resolution order:
    1. ``BIORUST_CLUSTALO_PATH`` env var (explicit path)
    2. Bundled platform binary under ``biorust/_bin/<platform>/clustalo``
    3. System PATH via ``shutil.which("clustalo")``
    """
    # 1. Env var override
    env = os.environ.get("BIORUST_CLUSTALO_PATH")
    if env is not None:
        p = Path(env)
        if not p.exists():
            raise FileNotFoundError(
                f"BIORUST_CLUSTALO_PATH points to '{p}' which does not exist"
            )
        if sys.platform != "win32" and not os.access(p, os.X_OK):
            raise PermissionError(
                f"BIORUST_CLUSTALO_PATH points to '{p}' which is not executable"
            )
        return p

    # 2. Bundled binary
    plat = None
    if sys.platform == "darwin" and platform.machine() in ("arm64", "aarch64"):
        plat = "darwin-arm64"
    elif sys.platform.startswith("linux") and platform.machine() == "x86_64":
        plat = "linux-x86_64"

    if plat is not None:
        ref = importlib.resources.files("biorust") / "_bin" / plat / "clustalo"
        path = Path(str(ref))
        if path.exists():
            os.chmod(path, 0o755)
            return path

    # 3. System PATH
    which = shutil.which("clustalo")
    if which is not None:
        return Path(which)

    raise RuntimeError(
        f"Could not find clustalo binary. Searched:\n"
        f"  1. BIORUST_CLUSTALO_PATH env var (not set)\n"
        f"  2. Bundled binary for {sys.platform}/{platform.machine()} (not found)\n"
        f"  3. System PATH via shutil.which (not found)\n"
        f"Install Clustal Omega or set BIORUST_CLUSTALO_PATH."
    )


def _normalize_records(records):
    """Convert input to ``(pairs, seq_type)`` where *seq_type* is ``"dna"`` or ``"protein"``."""
    from biorust import (
        DNARecord,
        DNARecordBatch,
        ProteinRecord,
        ProteinRecordBatch,
    )

    if isinstance(records, (str, bytes)):
        raise TypeError(
            "msa_clustalo() expects a RecordBatch or list[Record], not str or bytes"
        )

    if isinstance(records, DNARecordBatch):
        pairs = [(records[i].id, str(records[i].seq)) for i in range(len(records))]
        return pairs, "dna"

    if isinstance(records, ProteinRecordBatch):
        pairs = [(records[i].id, str(records[i].seq)) for i in range(len(records))]
        return pairs, "protein"

    if isinstance(records, list):
        if not records:
            raise ValueError("MSA requires at least 2 sequences")
        first = records[0]
        if isinstance(first, DNARecord):
            out = []
            for r in records:
                if not isinstance(r, DNARecord):
                    raise TypeError(
                        f"expected DNARecord in list, got {type(r).__name__}"
                    )
                out.append((r.id, str(r.seq)))
            return out, "dna"
        if isinstance(first, ProteinRecord):
            out = []
            for r in records:
                if not isinstance(r, ProteinRecord):
                    raise TypeError(
                        f"expected ProteinRecord in list, got {type(r).__name__}"
                    )
                out.append((r.id, str(r.seq)))
            return out, "protein"
        raise TypeError(
            f"expected DNARecord or ProteinRecord in list, got {type(first).__name__}"
        )

    raise TypeError(
        f"msa_clustalo() expects RecordBatch or list[Record], "
        f"got {type(records).__name__}"
    )


def _parse_fasta(text: str) -> list[tuple[str, str]]:
    """Parse simple FASTA text into (id, seq) pairs."""
    records = []
    current_id = None
    parts: list[str] = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                records.append((current_id, "".join(parts)))
            current_id = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)

    if current_id is not None:
        records.append((current_id, "".join(parts)))

    return records


def _validate_extra_args(extra_args):
    if isinstance(extra_args, (str, bytes)):
        raise TypeError("extra_args must be a list of strings, not str or bytes")
    if not isinstance(extra_args, (list, tuple)):
        raise TypeError(
            f"extra_args must be a list of strings, got {type(extra_args).__name__}"
        )

    forbidden_flags = {"-i", "-o", "--force"}
    forbidden_prefixes = (
        "--outfmt",
        "--seqtype",
        "--threads",
        "--in",
        "--out",
        "--infile",
        "--outfile",
    )

    for arg in extra_args:
        if not isinstance(arg, str):
            raise TypeError(
                f"extra_args must be a list of strings, got {type(arg).__name__}"
            )
        if arg in forbidden_flags or arg.startswith(forbidden_prefixes):
            raise ValueError(
                f"extra_args cannot override required clustalo args: {arg!r}"
            )


def msa_clustalo(records, *, threads=None, extra_args=None):
    """Run Clustal Omega multiple sequence alignment.

    Parameters
    ----------
    records : DNARecordBatch | list[DNARecord] | ProteinRecordBatch | list[ProteinRecord]
        Input sequences to align.
    threads : int, optional
        Number of threads for clustalo.
    extra_args : list[str], optional
        Additional command-line arguments for clustalo.

    Returns
    -------
    AlignmentDNA | AlignmentProtein
        Aligned (gapped) sequences.
    """
    from biorust import (
        AlignmentDNA,
        AlignmentProtein,
        GappedDNA,
        GappedProtein,
    )

    pairs, seq_type = _normalize_records(records)

    if len(pairs) < 2:
        raise ValueError("MSA requires at least 2 sequences")

    clustalo = find_clustalo()

    with tempfile.TemporaryDirectory() as tmpdir:
        in_path = os.path.join(tmpdir, "input.fasta")
        out_path = os.path.join(tmpdir, "output.fasta")

        with open(in_path, "w") as f:
            for seq_id, seq in pairs:
                f.write(f">{seq_id}\n{seq}\n")

        cmd = [
            str(clustalo),
            "-i",
            in_path,
            "-o",
            out_path,
            "--outfmt=fasta",
            "--force",
        ]

        if seq_type == "protein":
            cmd.append("--seqtype=Protein")

        if threads is not None:
            cmd.append(f"--threads={threads}")

        if extra_args:
            _validate_extra_args(extra_args)
            cmd.extend(extra_args)

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            stderr = result.stderr[:8192] if result.stderr else "(no stderr)"
            raise RuntimeError(
                f"clustalo failed (exit code {result.returncode})\n"
                f"command: {' '.join(cmd)}\n"
                f"stderr: {stderr}"
            )

        with open(out_path) as f:
            output_text = f.read()

    aligned = _parse_fasta(output_text)

    if not aligned:
        raise RuntimeError("clustalo produced no output sequences")

    width = len(aligned[0][1])
    for seq_id, seq in aligned:
        if len(seq) != width:
            raise RuntimeError(
                f"clustalo output has inconsistent lengths: "
                f"'{seq_id}' has {len(seq)}, expected {width}"
            )

    if seq_type == "protein":
        tuples = [(seq_id, GappedProtein(seq)) for seq_id, seq in aligned]
        return AlignmentProtein(tuples)

    tuples = [(seq_id, GappedDNA(seq)) for seq_id, seq in aligned]
    return AlignmentDNA(tuples)
