import pandas as pd

def make_incl_var_format(df):
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

    #The specified file has 5 columns with a header: RSID, position, chromosome, A_allele and B_allele                
    # the A and B allele thing needs to be checked by users
    # for finngen, they are swapped (A1, then A2) relative to 'rsids' (A2, then A1) just for fun I guess
    #head ~/meta_test.out
    #index RSID position chromosome A_allele B_allele A_allele_freq B_allele_freq
    #1 chr19_60842_A_G 60842 19 G A 0.0243708609 0.9756291391

    df = df[["id", "pos", "chrom", "A1", "A2"]]
    df.columns = ['RSID', 'position', 'chromosome', 'A_allele', 'B_allele']

    return df
    

sumstats_file = snakemake.input.snp_list
out = snakemake.output.out
print("Reading summary statistics", flush=True)

# Assumes tab-delimited input. Change sep if needed.
sumstats_chr = pd.read_csv(sumstats_file, sep="\t")

dat = sumstats_chr.copy()

dat = make_incl_var_format(dat)

dat.to_csv(out, sep=" ", index=False)
