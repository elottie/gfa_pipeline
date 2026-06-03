
## LD pruning options
l2_dir = config["analysis"]["R"]["l2_dir"]
ld_strings = expand("r2{r2}_kb{kb}_{p}",
                    r2 = config["analysis"]["ldprune"]["r2_thresh"],
                    kb = config["analysis"]["ldprune"]["clump_kb"],
                    p = config["analysis"]["ldprune"]["ld_prioritization"])
                    
## R options
if "pt" in config["analysis"]["R"]["type"]:
    R_strings = expand("pt{pt}",
                    pt = config["analysis"]["R"]["pthresh"])
else:
    R_strings = []

if "ldsc" in config["analysis"]["R"]["type"]:
    R_strings.append("ldsc")

if "none" in config["analysis"]["R"]["type"]:
    R_strings.append("none")

cor_cutoff = config["analysis"]["R"]["cor_cutoff"]
 
ldsc_mem_lim_gb = config["analysis"]["R"]["ldsc_mem_lim_gb"]    
    

# This produces one data frame per chromosome with columns for snp info
# and columns <study>.z, <study>.ss for z-score and sample size of each snp
## rule for getting list of raw data files
#def raw_data_input(wcs):
#    global prefix_dict
#    mycsv = prefix_dict[wcs.prefix]
#    ss = pd.read_csv(mycsv, na_filter=False)
#    return ss['raw_data_path']

def info_input(wcs):
    global prefix_dict
    return prefix_dict[wcs.prefix]

import re

def uncorr_info_input(wcs):
    return re.sub(r"\.csv$", "_uncorr_traits.csv", info_input(wcs))

# do we need 'raw data input' in this and the following rule?  removed from this rule
rule sample_size_bounds:
    input: gwas_info = uncorr_info_input
    output: out =  data_dir + "snp_lists/" + "{prefix}_sample_size_table.tsv"
    params: sample_size_tol = sstol_max    
    shell:  'bash bash/2_get_ss_bounds.sh {input.gwas_info} {params.sample_size_tol} {output.out}'

  
rule gather_snps:
    input: gwas_info = uncorr_info_input, 
           sample_size_file = data_dir + "snp_lists/" + "{prefix}_sample_size_table.tsv"
    output: out =  data_dir + "snp_lists/" + "{prefix}_snps_chr{chrom}.tsv" # output is two column file with list of rsids and minimum p-value/max z-score
    params: af_thresh = af_min,
            is_mvmr = is_mvmr
    wildcard_constraints: chrom = r"\d+"
    resources: mem_mb = 3000 # could adjust resources
            # add is_mvmr to script at some point
    shell: 'bash bash/2_gather_snps.sh {wildcards.chrom} {input.gwas_info} {input.sample_size_file} {params.af_thresh} {output.out}' 


# LD prune with plink
pthresh = 1 # jean change later or remove
rule ld_prune_plink:
    input: snp_list = data_dir + "snp_lists/" + "{prefix}_snps_chr{chrom}.tsv"
    output: out = data_dir + "snp_lists/" + "{prefix}_pruned_snps_r2{r2}_kb{kb}_{p}.{chrom}.tsv"
    params: ref_path = config["analysis"]["ldprune"]["ref_path"],
            pthresh = pthresh
    wildcard_constraints: chrom = r"\d+"
    resources: mem_mb = 5000 # could adjust resources
    script: 'R/3_ld_prune_chrom.R' # to update

# eventually needs diff options for non-GFA, ex. "beta" for beta and se for MRs
rule make_nice_data:
    input: gwas_info = uncorr_info_input,
           pruned_snp_list = expand(data_dir + "snp_lists/" + "{{prefix}}_pruned_snps_{{ldstring}}.{chrom}.tsv", chrom = range(1, 23))
    params: usage = "gfa"  # would be MR for those which want beta & se
    output: out = data_dir + "{prefix}_ldpruned_{ldstring}_nice_data_for_gfa.RData"
    script: "R/4_make_nice_data.R"


## Estimate R

# For p-value threshold and ldsc_quick methods, we can compute R
# without ever reading in all of the data.
# For ldsc method, we need to run ldsc for each pair of traits first.

### For strip division

checkpoint make_ldsc_strip_list:
    input: gwas_info = info_input
    output: out = data_dir + "{prefix}_ldsc_strip_list.json"
    params: mem_limit = ldsc_mem_lim_gb
    script: "R/1_R_make_ldsc_strip_list.R"


####p-value threshold method

#rule pt_R:
#  input: Z = expand(data_dir + "{{prefix}}_zmat.ldpruned_r2{{r2}}_kb{{kb}}_{{p}}.{chrom}.RDS", chrom = range(1, 23))
#  output: out = data_dir + "{prefix}_R_estimate.ldpruned_r2{r2}_kb{kb}_{p}.R_pt{pt}.RDS"
#  params: cond_num = cond_num
#  wildcard_constraints: pt = r"[\d.]+"
#  script: "R/3_R_pthresh.R"


### None
#rule none_R:
#    input: gwas_info = info_input
#    output: out = data_dir + "{prefix}_R_estimate.R_none.RDS"
#    script: 'R/3_R_none.R'


#rule R_ldsc_full:
#    input: Z = expand(data_dir + "{{prefix}}_zmat.{chrom}.RDS", chrom = range(1, 23)),
#           gwas_info = info_input,
#           m = expand(l2_dir + "{chrom}.l2.M_5_50", chrom = range(1, 23)),
#           l2 = expand(l2_dir + "{chrom}.l2.ldscore.gz", chrom = range(1, 23))
#    output: out = data_dir + "{prefix}_R_estimate.R_ldsc.RDS"
#    params: cond_num = cond_num
#    wildcard_constraints: pt = r"[\d.]+"
#    script: "R/3_R_ldsc_all.R"

# we need to ensure strip numbers passed in are from 1:length(strip_list-1).  cuz it will handle the last one (length(strip_list)) interally for free
rule R_ldsc_strip:
    input: gwas_info = info_input, 
           strip_list =  data_dir + "{prefix}_ldsc_strip_list.json",
           m = expand(l2_dir + "{chrom}.l2.M_5_50", chrom = range(1, 23)),
           l2 = expand(l2_dir + "{chrom}.l2.ldscore.gz", chrom = range(1, 23))
    output: out = data_dir + "{prefix}_R_estimate.R_ldsc.strip_{strip_num}.RDS"
    resources: mem_mb = ldsc_mem_lim_gb*1024 # could adjust resources
    script: "R/1_R_ldsc_strip.R"


# have to make this thing because we need to know number of strips from make_ldsc_strip_list (meaning it's a checkpoint)
# had to make strip list into json rather than rds so it could be read by python here
import json

def get_ldsc_strip_res(wcs):
    ck = checkpoints.make_ldsc_strip_list.get(prefix = wcs.prefix)
    json_file = ck.output.out

    with open(json_file) as f:
        ldsc_strips = json.load(f)

    nstrips = len(ldsc_strips)

    return expand(
        data_dir + "{prefix}_R_estimate.R_ldsc.strip_{strip_num}.RDS",
        prefix = wcs.prefix,
        strip_num = range(1, nstrips)
    )

rule R_ldsc_collect:
    input: gwas_info = info_input,
           ldsc_strip_res = get_ldsc_strip_res
    output: out = data_dir + "{prefix}_R_estimate.R_ldsc.RDS",
            uncorr_info = data_dir + "{prefix}_uncorr_traits.csv" 
    params: cor_cutoff = cor_cutoff,
            cond_num = cond_num
    script: "R/1_R_collect_ldsc_strips.R"

#import json
#with open(data_dir + f"{wcs.prefix}_ldsc_strip_list.json") as f:
#    ldsc_strips = json.load(f)
#nstrips = len(ldsc_strips)

#nstrips=2
#rule R_ldsc_collect:
#    input: ldsc_strip_res = expand(data_dir + "{{prefix}}_R_estimate.R_ldsc.{strip_num}.RDS", strip_num = range(1, nstrips))
#    output: out = data_dir + "{prefix}_R_estimate.R_ldsc.RDS"
#    script: "R/collect_ldsc_strips.R"
 






