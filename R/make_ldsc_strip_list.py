import math
import pandas as pd

# necessary not to use ceiling so we don't just get one block with <2 traits
def assign_traits(traits: list[str], nblocks: int) -> list[list[str]]:
    n = len(traits)
    base_num_traits = n // nblocks
    bigger_sets = n % nblocks

    out: list[list[str]] = []
    idx = 0
    for i in range(nblocks):
        set_size = base_num_traits + (1 if i < bigger_sets else 0)
        out.append(traits[idx:idx + set_size])
        idx += set_size
    return out

def make_trait_sets(
    csv_path: str,
    name_col: str = "name",
    memory_limit_gb: int = 150,
    trait_memory_gb: int = 5,
    blocks_at_once: int = 2,
    min_traits_per_block: int = 2,
) -> dict[str, list[str]]:

    traits = pd.read_csv(csv_path)[name_col].dropna().astype(str).tolist()
    ntraits = len(traits)

    if ntraits < min_traits_per_block:
        raise ValueError(f"Need at least {min_traits_per_block} traits; got {ntraits}.")

    max_traits_per_block = math.floor(memory_limit_gb / (blocks_at_once * trait_memory_gb))
    if max_traits_per_block < min_traits_per_block:
        raise ValueError(
            f"max_traits_per_block={max_traits_per_block} < min_traits_per_block={min_traits_per_block} "
            f"given memory_limit_gb={memory_limit_gb}, blocks_at_once={blocks_at_once}, trait_memory_gb={trait_memory_gb}."
        )

    # Feasible nblocks range
    nblocks_min = math.ceil(ntraits / max_traits_per_block)   # enough blocks to keep size <= max
    nblocks_max = ntraits // min_traits_per_block             # not too many blocks so size >= min

    if nblocks_min > nblocks_max:
        raise ValueError(
            "Impossible to satisfy block size constraints: "
            f"ntraits={ntraits}, min={min_traits_per_block}, max={max_traits_per_block}."
        )

    # Fewest blocks (goal)
    nblocks = nblocks_min

    strip_list = assign_traits(traits, nblocks)

    # Make absolutely sure no invalid blocks
    sizes = [len(v) for v in strip_list]
    if min(sizes) < min_traits_per_block or max(sizes) > max_traits_per_block:
        raise RuntimeError(f"Internal error: block sizes {sizes} violate [{min_traits_per_block}, {max_traits_per_block}].")

    return strip_list


if __name__ == "__main__":
    import json

    sets = make_trait_sets(csv_path=snakemake.input[0], memory_limit_gb=snakemake.params[0])
    with open(snakemake.output[0], "w") as f:
        json.dump(sets, f)

