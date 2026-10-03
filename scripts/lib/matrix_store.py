"""Unitig presence store and per-model pattern matrices.

Two levels, both written as a stream and never held whole in memory:

* The **organism store** holds every unitig that unitig-caller found in the
  organism's genomes, as one row of packed bits per unitig (bit j = genome j):

      genomes.csv          Genome ID, in bit order
      sequences.txt        line i = sequence of unitig i
      presence.bin         uint8, shape (n_unitigs, ceil(n_genomes / 8))
      store_summary.json   shapes, filter counts, source checksum

* A **model matrix** covers the genomes of one organism–antibiotic pair. Unitigs
  outside the frequency filter (present in fewer than ``min_support`` genomes or
  absent from fewer than ``min_support``) are dropped, and unitigs with identical
  presence/absence over the model's genomes are collapsed into one *pattern*:

      genomes.csv          Genome ID, label, lineage — row order of X
      X.npy                uint8, shape (n_genomes, ceil(n_patterns / 8)); bit p of
                           row g = pattern p present in genome g
      patterns.bin         uint8, shape (n_patterns, ceil(n_genomes / 8)); the same
                           matrix pattern by pattern
      patterns.csv         pattern_id, n_present, n_members
      members.npy          int64, shape (n_kept, 2): store unitig index, pattern id
      matrix_summary.json  shapes, filter counts, checksums

Packed bits use numpy's big-endian bit order within a byte (``np.packbits``); the
padding bits of the last byte are zero.
"""
from __future__ import annotations

import hashlib
import json
from collections.abc import Iterator
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

_ZERO, _ONE, _TAB = ord("0"), ord("1"), ord("\t")


def _nbytes(n_bits: int) -> int:
    return (n_bits + 7) // 8


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def _iter_lines(path: Path, block_bytes: int) -> Iterator[list[bytes]]:
    """Complete lines of a text file, in blocks of roughly ``block_bytes``."""
    with open(path, "rb") as f:
        rest = b""
        while True:
            chunk = f.read(block_bytes)
            if not chunk:
                break
            chunk = rest + chunk
            cut = chunk.rfind(b"\n")
            if cut < 0:
                rest = chunk
                continue
            rest = chunk[cut + 1:]
            yield chunk[:cut].split(b"\n")
        if rest.strip():
            yield [rest]


# ---------------------------------------------------------------------------
# organism store
# ---------------------------------------------------------------------------
def build_store(rtab: str | Path, out_dir: str | Path, *, min_support: int,
                block_bytes: int = 1 << 28) -> dict:
    """Parse unitig-caller's Rtab into the organism store.

    A unitig present in fewer than ``min_support`` genomes, or absent from fewer
    than ``min_support``, can pass the model filter in no subset of the genomes,
    so it is dropped here already. Lines are parsed a block at a time: the
    presence part of every line has the same width, so a whole block becomes one
    numpy array and no value is handled in Python.
    """
    rtab, out_dir = Path(rtab), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(rtab, "rb") as f:
        header = f.readline().rstrip(b"\r\n").split(b"\t")
    genomes = [g.decode() for g in header[1:]]
    n = len(genomes)
    if n == 0 or len(set(genomes)) != n:
        raise ValueError(f"{rtab}: the header must name each genome once")
    width = 2 * n - 1                        # "0\t1\t...\t0"
    pd.DataFrame({"Genome ID": genomes}).to_csv(out_dir / "genomes.csv", index=False)

    n_seen = n_kept = 0
    with open(out_dir / "presence.bin", "wb") as pres, \
            open(out_dir / "sequences.txt", "wb") as seqs:
        first = True
        for lines in _iter_lines(rtab, block_bytes):
            if first:                         # the header line
                lines, first = lines[1:], False
            lines = [ln.rstrip(b"\r") for ln in lines if ln.strip()]
            if not lines:
                continue
            tabs = [ln.find(b"\t") for ln in lines]
            body = [ln[t + 1:] for ln, t in zip(lines, tabs, strict=True)]
            if any(t < 1 or len(b) != width for t, b in zip(tabs, body, strict=True)):
                bad = next(i for i, (t, b) in enumerate(zip(tabs, body, strict=True))
                           if t < 1 or len(b) != width)
                raise ValueError(f"{rtab}: malformed line ({lines[bad][:60]!r}...): "
                                 f"expected {n} tab-separated 0/1 values")
            full = np.frombuffer(b"".join(body), dtype=np.uint8).reshape(len(lines), width)
            vals = full[:, ::2]
            if not (np.isin(vals, (_ZERO, _ONE)).all() and (full[:, 1::2] == _TAB).all()):
                raise ValueError(f"{rtab}: presence values must be tab-separated 0 or 1")
            present = vals == _ONE
            counts = present.sum(axis=1)
            keep = (counts >= min_support) & (n - counts >= min_support)
            n_seen += len(lines)
            if keep.any():
                pres.write(np.packbits(present[keep], axis=1).tobytes())
                seqs.write(b"".join(lines[i][:tabs[i]] + b"\n" for i in np.flatnonzero(keep)))
                n_kept += int(keep.sum())

    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_rtab": str(rtab), "source_sha256": sha256_file(rtab),
        "n_genomes": n, "n_unitigs_seen": n_seen, "n_unitigs": n_kept,
        "min_support": min_support,
        "presence_shape": [n_kept, _nbytes(n)],
        "tools": (json.loads((rtab.parent / "versions.json").read_text())
                  if (rtab.parent / "versions.json").exists() else {}),
    }
    (out_dir / "store_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


class Store:
    """Read access to an organism store."""

    def __init__(self, store_dir: str | Path):
        self.dir = Path(store_dir)
        self.summary = json.loads((self.dir / "store_summary.json").read_text())
        self.genomes = pd.read_csv(self.dir / "genomes.csv", dtype=str)["Genome ID"].tolist()
        shape = tuple(self.summary["presence_shape"])
        self.presence = (np.memmap(self.dir / "presence.bin", dtype=np.uint8, mode="r",
                                   shape=shape) if shape[0] else np.zeros(shape, np.uint8))

    @property
    def n_unitigs(self) -> int:
        return int(self.summary["n_unitigs"])

    def sequences(self, indices) -> list[str]:
        """Sequences of the given unitig indices, in the order given."""
        want = {int(i): k for k, i in enumerate(indices)}
        out: list[str | None] = [None] * len(want)
        with open(self.dir / "sequences.txt", encoding="ascii") as f:
            for i, line in enumerate(f):
                if i in want:
                    out[want[i]] = line.rstrip("\n")
        return [s for s in out if s is not None]


# ---------------------------------------------------------------------------
# model matrix
# ---------------------------------------------------------------------------
def build_model_matrix(store_dir: str | Path, genomes: pd.DataFrame, out_dir: str | Path, *,
                       min_support: int, block_unitigs: int = 1 << 16) -> dict:
    """Frequency filter and pattern collapse for one model.

    ``genomes`` has columns ``Genome ID``, ``label`` and ``lineage``; its order is
    the row order of X. Every genome must be in the store.
    """
    store, out_dir = Store(store_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    genomes = genomes[["Genome ID", "label", "lineage"]].reset_index(drop=True)
    pos = {g: i for i, g in enumerate(store.genomes)}
    missing = [g for g in genomes["Genome ID"] if g not in pos]
    if missing:
        raise ValueError(f"{len(missing)} genome(s) are not in the store, e.g. {missing[:3]}")
    cols = np.array([pos[g] for g in genomes["Genome ID"]], dtype=np.int64)
    n = len(cols)
    genomes.to_csv(out_dir / "genomes.csv", index=False)

    pattern_id: dict[bytes, int] = {}
    n_present: list[int] = []
    n_members: list[int] = []
    member_rows: list[np.ndarray] = []
    n_filtered = 0
    with open(out_dir / "patterns.bin", "wb") as pat:
        for start in range(0, store.n_unitigs, block_unitigs):
            stop = min(start + block_unitigs, store.n_unitigs)
            bits = np.unpackbits(np.asarray(store.presence[start:stop]), axis=1)[:, cols]
            counts = bits.sum(axis=1, dtype=np.int64)
            keep = np.flatnonzero((counts >= min_support) & (n - counts >= min_support))
            n_filtered += (stop - start) - keep.size
            if not keep.size:
                continue
            packed = np.packbits(bits[keep], axis=1)
            ids = np.empty(keep.size, dtype=np.int64)
            for k, row in enumerate(packed):
                key = row.tobytes()
                pid = pattern_id.get(key)
                if pid is None:
                    pid = pattern_id[key] = len(n_present)
                    n_present.append(int(counts[keep[k]]))
                    n_members.append(0)
                    pat.write(key)
                n_members[pid] += 1
                ids[k] = pid
            member_rows.append(np.column_stack([keep + start, ids]))
    n_patterns = len(n_present)
    del pattern_id
    if n_patterns == 0:
        raise ValueError(f"no unitig passes the frequency filter (min_support={min_support}, "
                         f"{n} genomes)")
    members = (np.concatenate(member_rows) if member_rows
               else np.empty((0, 2), dtype=np.int64))
    np.save(out_dir / "members.npy", members)
    pd.DataFrame({"pattern_id": np.arange(n_patterns), "n_present": n_present,
                  "n_members": n_members}).to_csv(out_dir / "patterns.csv", index=False)

    # genome-major copy: transpose 8·k patterns at a time
    x = np.lib.format.open_memmap(out_dir / "X.npy", mode="w+", dtype=np.uint8,
                                  shape=(n, _nbytes(n_patterns)))
    if n_patterns:
        pmap = np.memmap(out_dir / "patterns.bin", dtype=np.uint8, mode="r",
                         shape=(n_patterns, _nbytes(n)))
        step = 8 * 4096
        for p0 in range(0, n_patterns, step):
            p1 = min(p0 + step, n_patterns)
            block = np.unpackbits(np.asarray(pmap[p0:p1]), axis=1)[:, :n]   # (p, n)
            x[:, p0 // 8:_nbytes(p1)] = np.packbits(block.T, axis=1)
        del pmap
    x.flush()
    del x

    labels = genomes["label"].astype(int)
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "store_dir": str(store_dir), "store_source_sha256": store.summary["source_sha256"],
        "min_support": min_support,
        "n_genomes": n, "n_resistant": int((labels == 1).sum()),
        "n_susceptible": int((labels == 0).sum()),
        "n_unitigs_store": store.n_unitigs, "n_unitigs_filtered_out": int(n_filtered),
        "n_unitigs_kept": int(members.shape[0]), "n_patterns": n_patterns,
        "x_shape": [n, _nbytes(n_patterns)], "patterns_shape": [n_patterns, _nbytes(n)],
        "x_sha256": sha256_file(out_dir / "X.npy"),
    }
    (out_dir / "matrix_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


class ModelMatrix:
    """Read access to one model matrix."""

    def __init__(self, model_dir: str | Path):
        self.dir = Path(model_dir)
        self.summary = json.loads((self.dir / "matrix_summary.json").read_text())
        self.genomes = pd.read_csv(self.dir / "genomes.csv", dtype={"Genome ID": str})
        self.n_patterns = int(self.summary["n_patterns"])
        self.x = np.load(self.dir / "X.npy", mmap_mode="r")

    @property
    def n_genomes(self) -> int:
        return len(self.genomes)

    @property
    def labels(self) -> np.ndarray:
        return self.genomes["label"].to_numpy(dtype=np.int8)

    def rows(self, row_indices) -> np.ndarray:
        """Dense uint8 0/1 block (len(row_indices) × n_patterns), rows in the order given."""
        idx = np.asarray(row_indices, dtype=np.int64)
        return np.unpackbits(np.asarray(self.x[idx]), axis=1)[:, :self.n_patterns]

    def row_batches(self, row_indices, batch_size: int) -> Iterator[np.ndarray]:
        idx = np.asarray(row_indices, dtype=np.int64)
        for s in range(0, idx.size, batch_size):
            yield self.rows(idx[s:s + batch_size])

    def pattern(self, pattern_id: int) -> np.ndarray:
        """0/1 column of one pattern over all genomes."""
        pmap = np.memmap(self.dir / "patterns.bin", dtype=np.uint8, mode="r",
                         shape=tuple(self.summary["patterns_shape"]))
        return np.unpackbits(np.asarray(pmap[pattern_id]))[:self.n_genomes]

    def present_in(self, genome_mask, block: int = 1 << 15) -> np.ndarray:
        """Number of genomes in ``genome_mask`` (one bool per genome) carrying each pattern."""
        mask = np.asarray(genome_mask, dtype=bool)
        if mask.shape != (self.n_genomes,):
            raise ValueError("genome_mask needs one value per genome")
        pmap = np.memmap(self.dir / "patterns.bin", dtype=np.uint8, mode="r",
                         shape=tuple(self.summary["patterns_shape"]))
        out = np.empty(self.n_patterns, dtype=np.int64)
        for s in range(0, self.n_patterns, block):
            bits = np.unpackbits(np.asarray(pmap[s:s + block]), axis=1)[:, :self.n_genomes]
            out[s:s + block] = bits[:, mask].sum(axis=1, dtype=np.int64)
        return out

    def columns(self, pattern_ids) -> np.ndarray:
        """Dense uint8 0/1 block (n_genomes × len(pattern_ids)) of the given patterns."""
        ids = np.asarray(pattern_ids, dtype=np.int64)
        if ids.size and (ids.min() < 0 or ids.max() >= self.n_patterns):
            raise IndexError(f"pattern ids must lie in [0, {self.n_patterns})")
        pmap = np.memmap(self.dir / "patterns.bin", dtype=np.uint8, mode="r",
                         shape=tuple(self.summary["patterns_shape"]))
        bits = np.unpackbits(np.asarray(pmap[ids]), axis=1)[:, :self.n_genomes]
        return np.ascontiguousarray(bits.T)

    def members(self) -> pd.DataFrame:
        """Store unitig index and pattern id of every kept unitig."""
        m = np.load(self.dir / "members.npy")
        return pd.DataFrame({"unitig_index": m[:, 0], "pattern_id": m[:, 1]})
