from make_ldsc_strip_list import make_trait_sets  # adjust path/name as needed
import json

#def make_trait_sets(
#    csv_path: str,
#    name_col: str = "name",
#    memory_limit_gb: int = 150,
#    trait_memory_gb: int = 5,
#    blocks_at_once: int = 2,
#    min_traits_per_block: int = 2,
#) -> dict[str, list[str]]:

#"../C100001554_And_Friends_3Metabolites.csv"
#"../5e5Sig_Herit_Mets_8ForLDSCStrip.csv"
sets = make_trait_sets(csv_path="../5e5Sig_Herit_Mets_8ForLDSCStrip.csv",memory_limit_gb=5,trait_memory_gb=1)
print(sets)

with open("ldsc_trait_sets.json", "w") as f:
    json.dump(sets, f)
