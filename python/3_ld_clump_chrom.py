# tell ppl to install ldstore w/ pip from christian benner's website into the conda env used for gfa

import numpy as np
import pandas as pd
from bisect import bisect_left, bisect_right
from ldstore.bcor import bcor


def normalize_chr(x):
    """
    Normalize chromosome labels:
      1, '1', '01', 'chr1', 'CHR1' -> '1'
    """
    x = str(x).replace("chr", "").replace("CHR", "")
    try:
        return str(int(x))
    except ValueError:
        return x


def ld_clump_chr(
    sumstats_chr,
    bcor_file,
    r2_threshold,
    distance_kb,
):
    """
    P-value-prioritized LD clumping for one chromosome using an LDstore .bcor file.

    Algorithm:
      1. Sort SNPs by p-value.
      2. Take the lowest-p SNP not already removed.
      3. Keep it as a lead SNP.
      4. Find nearby SNPs within +/- distance_kb.
      5. Use .bcor LD to remove SNPs with r^2 > r2_threshold.
      6. Repeat until all SNPs are kept or removed.
    """

    distance_bp = int(distance_kb * 1000)

    print(f"Opening BCOR file: {bcor_file}", flush=True)
    b = bcor(bcor_file)

    print("Reading BCOR metadata", flush=True)
    meta = b.getMeta().copy()

    # The row index in metadata corresponds to the row/column index in the .bcor LD matrix.
    meta["bcor_index"] = meta.index.astype(int)

    # Standardize BCOR metadata column names.
    meta = meta.rename(
        columns={
            "rsid": "snp",
            "position": "pos",
            "chromosome": "chr",
        }
    )

    required_meta_cols = {"snp", "pos", "chr", "bcor_index"}
    missing_meta = required_meta_cols - set(meta.columns)
    if missing_meta:
        raise ValueError(f"BCOR metadata missing required columns: {missing_meta}")

    meta["chr"] = meta["chr"].apply(normalize_chr)
    meta["pos"] = pd.to_numeric(meta["pos"], errors="coerce")
    meta = meta.dropna(subset=["snp", "pos", "chr"])
    meta["pos"] = meta["pos"].astype(int)

    print(f"BCOR contains {len(meta)} SNPs", flush=True)

    # Merge sumstats with BCOR metadata by SNP ID.
    dat = sumstats_chr.merge(
        meta[["snp", "chr", "pos", "bcor_index"]],
        on="snp",
        how="inner",
        suffixes=("", "_bcor"),
    )

    print(f"Matched {len(dat)} summary-stat SNPs to BCOR SNPs", flush=True)

    if dat.empty:
        return pd.DataFrame(), pd.DataFrame()

    # don't need to turn z-scores to pvals, and helps us avoid precision issues
    dat["max_abs_z"] = pd.to_numeric(dat["max_abs_z"], errors="coerce")
    dat = dat.dropna(subset=["max_abs_z"])

    if dat.empty:
        return pd.DataFrame(), pd.DataFrame()

    # Sort SNPs by p-value ascending, i.e. absolute Z-score descending.
    dat = dat.sort_values("max_abs_z", ascending=False).reset_index(drop=True)

    # A position-sorted version for fast window lookup.
    pos_sorted = dat.sort_values("pos").reset_index(drop=True)
    positions = pos_sorted["pos"].to_numpy()

    removed_snps = set()
    kept_records = []
    removed_records = []

    for _, lead in dat.iterrows():
        lead_snp = lead["snp"]

        if lead_snp in removed_snps:
            continue

        lead_pos = int(lead["pos"])
        lead_z = float(lead["max_abs_z"])
        lead_index = int(lead["bcor_index"])

        kept_records.append(
            {
                "lead_snp": lead_snp,
                "lead_pos": lead_pos,
                "lead_z": lead_z,
                "lead_bcor_index": lead_index,
            }
        )

        # Identify SNPs within +/- distance_bp.
        left = bisect_left(positions, lead_pos - distance_bp)
        right = bisect_right(positions, lead_pos + distance_bp)

        window = pos_sorted.iloc[left:right].copy()

        # Skip SNPs already removed by previous lead SNPs.
        window = window[~window["snp"].isin(removed_snps)]

        if window.empty:
            removed_snps.add(lead_snp)
            continue

        other_snps = window["snp"].to_numpy()
        other_pos = window["pos"].to_numpy()
        other_z = window["max_abs_z"].to_numpy()
        other_indices = window["bcor_index"].to_numpy(dtype=int)

        # readCorr([lead_index]) returns correlations involving the lead SNP.
        # Based on your example:
        #   myBcor.readCorr([29])[10, 0]
        # gives corr between SNP 10 and SNP 29.
        corr = np.asarray(b.readCorr([lead_index]))

        if corr.ndim == 1:
            r = corr[other_indices]
        else:
            r = corr[other_indices, 0]

        r2 = r ** 2
        high_ld = r2 > r2_threshold

        for snp, pos, z, idx, r2_val in zip(
            other_snps[high_ld],
            other_pos[high_ld],
            other_z[high_ld],
            other_indices[high_ld],
            r2[high_ld],
        ):
            removed_snps.add(snp)

            removed_records.append(
                {
                    "lead_snp": lead_snp,
                    "lead_pos": lead_pos,
                    "lead_z": lead_z,
                    "removed_snp": snp,
                    "removed_pos": int(pos),
                    "removed_z": float(z),
                    "removed_bcor_index": int(idx),
                    "r2": float(r2_val),
                }
            )

        # Make sure the lead SNP is marked as processed.
        removed_snps.add(lead_snp)

        if len(kept_records) % 1000 == 0:
            print(f"Kept {len(kept_records)} lead SNPs so far", flush=True)

    kept = pd.DataFrame(kept_records)
    removed = pd.DataFrame(removed_records)

    return kept, removed


print("Starting one-chromosome BCOR LD clumping", flush=True)

# ---------------------------------------------------------------------
# Inputs/outputs/params from Snakemake
# ---------------------------------------------------------------------

sumstats_file = snakemake.input.snp_list
bcor_file = snakemake.input.bcor_file

kept_out = snakemake.output.kept
removed_out = snakemake.output.removed
clumped_list_out = snakemake.output.clumped_snp_list

# Get chromosome from wildcard if available.
chrom = normalize_chr(snakemake.wildcards.chrom)

r2_threshold = float(snakemake.wildcards.r2)
distance_kb = float(snakemake.wildcards.kb)

print(f"Chromosome: {chrom}", flush=True)
print(f"Summary statistics: {sumstats_file}", flush=True)
print(f"BCOR file: {bcor_file}", flush=True)
print(f"r2 threshold: {r2_threshold}", flush=True)
print(f"distance kb: {distance_kb}", flush=True)

# ---------------------------------------------------------------------
# Read and clean summary statistics
# ---------------------------------------------------------------------

print("Reading summary statistics", flush=True)

# Assumes tab-delimited input. Change sep if needed.
sumstats = pd.read_csv(sumstats_file, sep="\t")

required_cols = {"snp", "max_abs_z"}
missing_cols = required_cols - set(sumstats.columns)

if missing_cols:
    raise ValueError(f"Summary statistics missing required columns: {missing_cols}")

sumstats = sumstats.copy()
sumstats["max_abs_z"] = pd.to_numeric(sumstats["max_abs_z"], errors="coerce")
sumstats = sumstats.dropna(subset=["snp", "max_abs_z"])

# Already filtered to this chromosome.

print(f"SNPs on chr{chrom} after cleaning: {len(sumstats)}", flush=True)

# ---------------------------------------------------------------------
# Process chromosome
# ---------------------------------------------------------------------

all_kept = []
all_removed = []

if sumstats.empty:
    print(f"No SNPs for chr{chrom}; writing empty outputs", flush=True)

    kept_empty = pd.DataFrame(
        columns=["chr", "lead_snp", "lead_pos", "lead_z", "lead_bcor_index"]
    )

    removed_empty = pd.DataFrame(
        columns=[
            "chr",
            "lead_snp",
            "lead_pos",
            "lead_z",
            "removed_snp",
            "removed_pos",
            "removed_z",
            "removed_bcor_index",
            "r2",
        ]
    )

    kept_empty.to_csv(kept_out, sep="\t", index=False)
    removed_empty.to_csv(removed_out, sep="\t", index=False)
    open(clumped_list_out, "w").close()

else:
    kept, removed = ld_clump_chr(
        sumstats_chr=sumstats,
        bcor_file=bcor_file,
        r2_threshold=r2_threshold,
        distance_kb=distance_kb,
    )

    if kept.empty:
        kept = pd.DataFrame(
            columns=["lead_snp", "lead_pos", "lead_z", "lead_bcor_index"]
        )

    if removed.empty:
        removed = pd.DataFrame(
            columns=[
                "lead_snp",
                "lead_pos",
                "lead_z",
                "removed_snp",
                "removed_pos",
                "removed_z",
                "removed_bcor_index",
                "r2",
            ]
        )

    kept.insert(0, "chr", chrom)
    removed.insert(0, "chr", chrom)

    print(f"Writing kept SNP table: {kept_out}", flush=True)
    kept.to_csv(kept_out, sep="\t", index=False)

    print(f"Writing removed SNP table: {removed_out}", flush=True)
    removed.to_csv(removed_out, sep="\t", index=False)

    print(f"Writing retained SNP list: {clumped_list_out}", flush=True)

    if not kept.empty:
        kept["lead_snp"].to_csv(clumped_list_out, index=False, header=False)
    else:
        open(clumped_list_out, "w").close()

print("One-chromosome BCOR LD clumping complete", flush=True)
