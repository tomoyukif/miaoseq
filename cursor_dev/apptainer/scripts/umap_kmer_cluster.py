#!/usr/bin/env python3
"""Coarse clustering of amplicon reads by 5-mer frequency + UMAP + HDBSCAN."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import hdbscan
import numpy as np
import umap
from sklearn.preprocessing import normalize


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    rid = None
    seq_parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if rid is not None:
                    records.append((rid, "".join(seq_parts).upper()))
                rid = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if rid is not None:
        records.append((rid, "".join(seq_parts).upper()))
    return records


def kmer_vector(seq: str, k: int = 5) -> Counter[str]:
    counts: Counter[str] = Counter()
    if len(seq) < k:
        return counts
    for i in range(len(seq) - k + 1):
        counts[seq[i : i + k]] += 1
    return counts


def build_matrix(records: list[tuple[str, str]], k: int = 5) -> tuple[list[str], np.ndarray]:
    vectors = [kmer_vector(seq, k=k) for _, seq in records]
    vocab = sorted({kmer for vec in vectors for kmer in vec})
    if not vocab:
        ids = [rid for rid, _ in records]
        return ids, np.zeros((len(ids), 1), dtype=float)

    mat = np.zeros((len(records), len(vocab)), dtype=float)
    index = {kmer: i for i, kmer in enumerate(vocab)}
    for row, vec in enumerate(vectors):
        for kmer, count in vec.items():
            mat[row, index[kmer]] = count
    mat = normalize(mat, norm="l1", axis=1)
    ids = [rid for rid, _ in records]
    return ids, mat


def cluster_reads(ids: list[str], matrix: np.ndarray, min_samples: int) -> list[tuple[str, int]]:
    if len(ids) == 0:
        return []
    if len(ids) == 1:
        return [(ids[0], 0)]

  # Small samples: one coarse group.
    if len(ids) < max(10, min_samples):
        return [(rid, 0) for rid in ids]

    n_neighbors = min(15, max(2, len(ids) - 1))
    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=0.0,
        metric="cosine",
        random_state=42,
    )
    embedding = reducer.fit_transform(matrix)
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=max(2, min_samples),
        min_samples=max(1, min_samples // 2),
        metric="euclidean",
    )
    labels = clusterer.fit_predict(embedding)
    return list(zip(ids, labels.astype(int)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta")
    parser.add_argument("output_tsv")
    parser.add_argument("--min-samples", type=int, default=50)
    parser.add_argument("--kmer", type=int, default=5)
    args = parser.parse_args()

    records = read_fasta(Path(args.fasta))
    ids, matrix = build_matrix(records, k=args.kmer)
    assignments = cluster_reads(ids, matrix, min_samples=args.min_samples)

    out = Path(args.output_tsv)
    with out.open("w") as handle:
        for rid, group in assignments:
            handle.write(f"{rid}\t{group}\n")


if __name__ == "__main__":
    main()
