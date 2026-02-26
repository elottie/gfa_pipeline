#!/bin/bash

# overall usage:
#source get_col.sh                # OR . get_col.sh

#parse_header "yourfile.tsv"      # Build the associative array col_indices

#snp_col=$(get_col "SNP")         # Retrieve the index for column "SNP"

#awk -v snp_col="$snp_col" ...

# Function to parse header and build associative array
parse_header() {
    local file="$1"

    if [[ ! -r "$file" ]]; then
        echo "Error: File '$file' does not exist or is not readable." >&2
        return 1
    fi

    if [[ "$file" == *.gz ]]; then
        header=$(awk 'NR==1 {print; exit}' <(gzip -cd "$file"))
    else
        header=$(awk 'NR==1 {print; exit}' "$file")
    fi

    if [[ -z "$header" ]]; then
        echo "Error: Header is empty or could not be read from '$file'." >&2
        return 1
    fi

    IFS=$'\t' read -r -a cols <<< "$header"
    declare -gA col_indices
    for i in "${!cols[@]}"; do
        colname=${cols[$i]}
        colname=$(echo "$colname" | sed 's/^[`"'\'' ]*//;s/[`"'\'' ]*$//')
        col_indices[$colname]=$((i+1))
    done
    return 0
}
# Usage: parse_header "myfile.tsv" (sets global col_indices)

get_col() {
    local name="$1"
    if [[ "$name" != "" && "$name" != "NA" ]]; then
        if [[ -v col_indices[$name] ]]; then
            echo "${col_indices[$name]}"
        else
            echo "Warning: column '$name' not found in header" >&2
            echo ""
        fi
    else
        echo ""
    fi
}

# Optional: main section for standalone use
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    file="$1"
    if [[ -z "$file" ]]; then
        echo "Usage: $0 <filename.tsv>" >&2
        exit 1
    fi
    parse_header "$file"
    echo "Column index for SNP:" $(get_col "SNP")
    echo "Column indices available:" "${!col_indices[@]}"
fi
