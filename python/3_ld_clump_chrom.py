import os
import subprocess
import tempfile
import numpy as np
import pandas as pd
import re
from bisect import bisect_left, bisect_right


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

def add_snp_id(df):
    """
    Add SNP ID column in format:
    chrCHROM_POS_ALT_REF

    Uses:
    - df["chrom"]
    - df["pos"]
    - df["A2"]  <- alt
    - df["A1"]  <- ref
    """

    df["id"] = (
        "chr" + df["chrom"].astype(str) + "_" +
        df["pos"].astype(str) + "_" +
        df["A2"].astype(str) + "_" +
        df["A1"].astype(str)
    )

    return df[["snp", "id", "pos", "max_abs_z"]]

def run_ldstore_range(ldstore_exec, bcor_file, start_bp, end_bp, incl_var_file):
    """
    Use LDstore CLI to extract LD for one genomic range.

    Equivalent shell command:
      ldstore --bcor FILE.bcor --incl-range START-END --table TEMP.tab

    Returns
    -------
    pandas.DataFrame
        Table produced by LDstore.
    """

    start_bp = max(0, int(start_bp))
    end_bp = int(end_bp)

    with tempfile.NamedTemporaryFile(
        suffix=".ldstore.tab",
        delete=False
    ) as tmp:
        tmp_file = tmp.name

    # could add ld-thold?
    # ARGS MUST BE MATCHING ORDER OF HELP AND IT WONT TELL YOU THAT
    cmd = [
        ldstore_exec,
        "--bcor", str(bcor_file),
        "--table", str(tmp_file),
        #"--incl-range", f"{start_bp}-{end_bp}",
        "--incl-variants", str(incl_var_file),
    ]


    with tempfile.NamedTemporaryFile(
        suffix=".ldstore.tab",
        delete=False
    ) as tmp:
        tmp_file_2 = tmp.name

    #cmd_2 = [
    #    ldstore_exec,
    #    "--bcor", str(bcor_file),
    #    "--incl-range", f"{start_bp}-{end_bp}",
    #    "--meta", str(tmp_file_2),
    #    #"--incl-variants", str(incl_var_file),
    #]

    try:
        print("Running LDstore:", " ".join(cmd), flush=True)

        result = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        #result_2 = subprocess.run(
        #    cmd_2,
        #    check=False,
        #    stdout=subprocess.PIPE,
        #    stderr=subprocess.PIPE,
        #    text=True,
        #)

        if result.stdout:
            print(result.stdout, flush=True)

        if result.stderr:
            print(result.stderr, flush=True)

        if not os.path.exists(tmp_file) or os.path.getsize(tmp_file) == 0:
            return pd.DataFrame()

#        import shlex

        #if result.returncode != 0:
        #    raise RuntimeError(
        #    f"STDOUT:\n{result.stdout}\n\n"
        #    f"STDERR:\n{result.stderr}"
        #)

 #       if result_2.returncode != 0:
 #           raise RuntimeError(
 #           f"LDstore failed with return code {result.returncode}\n\n"
 #           f"Command:\n{shlex.join(cmd)}\n\n"
 #           f"STDOUT:\n{result.stdout}\n\n"
 #           f"STDERR:\n{result.stderr}"
 #       )

        # LDstore table output is usually whitespace-delimited.
        # If your inspected file is tab-delimited only, sep="\t" is also fine.
        ld_table = pd.read_csv(tmp_file, sep=r"\s+", engine="python")
        
        print('Successfully fetched ld_table:')

        print(ld_table.head())

    finally:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)

    return ld_table

def check_ldstore_ids(
    ld_table,
    cols=("RSID1", "RSID2"),
    max_examples=10
):
    """
    Check that all IDs in ldstore rsid columns are of the form:

        chrCHROM_POS_ALT_REF

    Example:
        chr1_12345_A_G
    """
    missing_cols = [c for c in cols if c not in ld_table.columns]
    if missing_cols:
        raise ValueError(
            f"ldstore table is missing required column(s): {missing_cols}. "
            f"Available columns are: {list(ld_table.columns)}"
        )

    
    bad_row_mask = False
    messages = []

    SNP_ID_PATTERN = re.compile(
        r"^chr(?P<chrom>[1-9]|1[0-9]|2[0-2])_(?P<pos>[1-9][0-9]*)_(?P<alt>[ACGT]+)_(?P<ref>[ACGT]+)$"
    )

    for col in cols:
        values = ld_table[col].astype("string")
        format_ok = values.str.fullmatch(SNP_ID_PATTERN, na=False)

        bad_mask = ~format_ok

        if bad_mask.any():
            bad_row_mask = bad_row_mask | bad_mask
            messages.append(
                f"Column '{col}' has {bad_mask.sum()} invalid SNP IDs."
            )

    if isinstance(bad_row_mask, bool):
        return

    if bad_row_mask.any():
        examples = ld_table.loc[bad_row_mask, list(cols)].head(max_examples)

        raise ValueError(
            "Invalid ldstore SNP IDs detected.\n"
            "Expected all values in rsid1 and rsid2 to match:\n"
            "    chrCHROM_POS_ALT_REF\n\n"
            + "\n".join(messages)
            + "\n\nExample offending rows:\n"
            + examples.to_string(index=True)
        )
    else:
        print('No invalid ldstore SNP IDs detected')

def extract_lead_ld(ld_table, lead_id, r2_threshold):
    """
    Parse LDstore v1.1 --table output with columns:

    chromosome index1 RSID1 position1 index2 RSID2 position2 correlation n_samples

    Return rows for SNPs with r^2 > r2_threshold involving lead_id.
    """

    required = {
        "RSID1",
        "position1",
        "RSID2",
        "position2",
        "correlation",
    }

    missing = required - set(ld_table.columns)

    if missing:
        raise ValueError(
            f"LDstore table missing required columns: {missing}. "
            f"Observed columns: {list(ld_table.columns)}"
        )

    df = ld_table.copy()

    df["RSID1"] = df["RSID1"].astype(str)
    df["RSID2"] = df["RSID2"].astype(str)
    df["position1"] = pd.to_numeric(df["position1"], errors="coerce")
    df["position2"] = pd.to_numeric(df["position2"], errors="coerce")
    df["correlation"] = pd.to_numeric(df["correlation"], errors="coerce")

    df = df.dropna(subset=["RSID1", "RSID2", "position1", "position2", "correlation"])

    # Rows involving lead SNP.
    hit = df[(df["RSID1"] == str(lead_id)) | (df["RSID2"] == str(lead_id))].copy()

    if hit.empty:
        return pd.DataFrame(columns=["removed_id", "removed_pos", "r2"])

    hit["r2"] = hit["correlation"] ** 2
    hit = hit[hit["r2"] > r2_threshold].copy()

    if hit.empty:
        return pd.DataFrame(columns=["removed_id", "removed_pos", "r2"])

    # Determine the other SNP and its position.
    lead_is_1 = hit["RSID1"] == str(lead_id)

    hit["removed_id"] = np.where(
        lead_is_1,
        hit["RSID2"],
        hit["RSID1"],
    )

    hit["removed_pos"] = np.where(
        lead_is_1,
        hit["position2"],
        hit["position1"],
    )

    print('Successfully returning lead_ld results:  removed_id, removed_pos, r2')

    return hit[["removed_id", "removed_pos", "r2"]].drop_duplicates()

def ld_clump_chr(
    sumstats_chr,
    incl_var_file,
    bcor_file,
    ldstore_exec,
    r2_threshold,
    distance_kb,
):
    """
    Z-score-prioritized LD clumping for one chromosome using LDstore CLI.

    Algorithm:
      1. Sort SNPs by max_abs_z descending.
      2. Take the highest-Z SNP not already excluded.
      3. Keep it as a lead SNP.
      4. Find nearby SNPs within +/- distance_kb using positions in sumstats_chr.
      5. Use LDstore CLI to extract LD for lead_pos +/- distance_kb.
      6. Remove SNPs with r^2 > r2_threshold.
      7. Repeat until all SNPs are kept or excluded.

    Requirements
    ------------
    sumstats_chr must contain:
      - chrom
      - snp
      - pos
      - A2  <- alt
      - A1  <- ref
      - max_abs_z
    """

    distance_bp = int(distance_kb * 1000)

    print(f"Using LDstore executable: {ldstore_exec}", flush=True)
    print(f"Using BCOR file: {bcor_file}", flush=True)

    dat = sumstats_chr.copy()

    required_cols = {"chrom", "snp", "pos", "A2", "A1", "max_abs_z"}
    missing_cols = required_cols - set(dat.columns)

    if missing_cols:
        raise ValueError(
            f"sumstats_chr is missing required columns: {missing_cols}. "
            "With the LDstore CLI approach, the SNP list must include positions "
            "because Python can no longer read BCOR metadata directly."
        )


    dat = add_snp_id(dat)
    dat["pos"] = pd.to_numeric(dat["pos"], errors="coerce")
    dat["max_abs_z"] = pd.to_numeric(dat["max_abs_z"], errors="coerce")

    dat = dat.dropna(subset=["snp", "id", "pos", "max_abs_z"])
    dat["pos"] = dat["pos"].astype(int)

    if dat.empty:
        return pd.DataFrame(), pd.DataFrame()

    # No duplicate SNPs exist bc of preprocessing.
    dat = dat.sort_values("max_abs_z", ascending=False).reset_index(drop=True)

    print(f"After cleaning/deduplication: {len(dat)} SNPs", flush=True)

    # Position-sorted table for fast local-window lookup.
    pos_sorted = dat.sort_values("pos").reset_index(drop=True)
    positions = pos_sorted["pos"].to_numpy()

    excluded_snps = set()
    kept_records = []
    removed_records = []

    for _, lead in dat.iterrows():
        lead_snp = lead["snp"]

        if lead_snp in excluded_snps:
            continue

        lead_id = lead["id"]
        lead_pos = int(lead["pos"])
        lead_z = float(lead["max_abs_z"])

        kept_records.append(
            {
                "lead_snp": lead_snp,
                "lead_id": lead_id,
                "lead_pos": lead_pos,
                "lead_z": lead_z,
            }
        )

        # Identify SNPs within +/- distance_bp based on the SNP-list positions.
        left = bisect_left(positions, lead_pos - distance_bp)
        right = bisect_right(positions, lead_pos + distance_bp)

        window = pos_sorted.iloc[left:right].copy()

        # Skip SNPs already kept or removed by previous lead SNPs.
        window = window[~window["snp"].isin(excluded_snps)]

        if window.empty:
            excluded_snps.add(lead_snp)
            continue

        # Extract LD for this lead SNP's local range using LDstore CLI.
        ld_table = run_ldstore_range(
            ldstore_exec=ldstore_exec,
            bcor_file=bcor_file,
            start_bp=lead_pos - distance_bp,
            end_bp=lead_pos + distance_bp,
            incl_var_file=incl_var_file,
        )

        check_ldstore_ids(ld_table)

        # Get SNPs in LD with the lead.
        high_ld = extract_lead_ld(
            ld_table=ld_table,
            lead_id=lead_id,
            r2_threshold=r2_threshold,
        )

        # Restrict removals to SNPs in our current candidate window/list.
        window_info = window[["snp", "id", "max_abs_z"]].copy()

        high_ld = high_ld.merge(
            window_info,
            left_on="removed_id",
            right_on="id",
            how="inner",
        )

        # Always exclude the lead itself from future consideration.
        excluded_snps.add(lead_snp)

        for _, row in high_ld.iterrows():
            removed_snp = row["snp"]

            excluded_snps.add(removed_snp)

            # Do not record the lead SNP as removed.
            if removed_snp == lead_snp:
                continue

            removed_records.append(
                {
                    "lead_snp": lead_snp,
                    "lead_id": lead_id,
                    "lead_pos": lead_pos,
                    "lead_z": lead_z,
                    "removed_snp": removed_snp,
                    "removed_id": row["removed_id"],
                    "removed_pos": int(row["removed_pos"]),
                    "removed_z": float(row["max_abs_z"]),
                    "r2": float(row["r2"]),
                }
            )

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
incl_var_file = snakemake.input.incl_var_file
bcor_file = snakemake.input.bcor_file
ldstore_exec = snakemake.params.ldstore_exec

kept_out = snakemake.output.kept
removed_out = snakemake.output.removed
clumped_list_out = snakemake.output.clumped_snp_list

# Get chromosome from wildcard if available.
chrom = normalize_chr(snakemake.wildcards.chrom)

r2_threshold = float(snakemake.wildcards.r2)
distance_kb = float(snakemake.wildcards.kb)

print(f"Chromosome: {chrom}", flush=True)
print(f"Summary statistics: {sumstats_file}", flush=True)
print(f"Included var file for ldstore: {incl_var_file}", flush=True)
print(f"BCOR file: {bcor_file}", flush=True)
print(f"LDstore executable file: {ldstore_exec}", flush=True)
print(f"r2 threshold: {r2_threshold}", flush=True)
print(f"distance kb: {distance_kb}", flush=True)

# ---------------------------------------------------------------------
# Read and clean summary statistics
# ---------------------------------------------------------------------

print("Reading summary statistics", flush=True)

# Assumes tab-delimited input. Change sep if needed.
sumstats = pd.read_csv(sumstats_file, sep="\t")

# Already filtered to this chromosome.
print(f"SNPs on chr{chrom} after cleaning: {len(sumstats)}", flush=True)

# ---------------------------------------------------------------------
# Process chromosome
# ---------------------------------------------------------------------

if sumstats.empty:
    print(f"No SNPs for chr{chrom}; writing empty outputs", flush=True)

    kept_empty = pd.DataFrame(
        columns=["chr", "lead_snp", "lead_id", "lead_pos", "lead_z"]
    )

    removed_empty = pd.DataFrame(
        columns=[
            "chr",
            "lead_snp",
            "lead_id",
            "lead_pos",
            "lead_z",
            "removed_snp",
            "removed_id",
            "removed_pos",
            "removed_z",
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
        ldstore_exec=ldstore_exec,
        r2_threshold=r2_threshold,
        distance_kb=distance_kb,
        incl_var_file=incl_var_file,
    )

    if kept.empty:
        kept = pd.DataFrame(
            columns=["lead_snp", "lead_id", "lead_pos", "lead_z"]
        )

    if removed.empty:
        removed = pd.DataFrame(
            columns=[
                "lead_snp",
                "lead_id",
                "lead_pos",
                "lead_z",
                "removed_snp",
                "removed_id",
                "removed_pos",
                "removed_z",
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

stop('you told me to stop here')
