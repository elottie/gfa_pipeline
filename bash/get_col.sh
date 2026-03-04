#!/bin/bash

# purpose:
# a set of functions to get the column indices given the column names, since this is not easy in bash

# usage:
# source get_col.sh
# parse_header "$file" "$delimiter"      # build the associative array col_indices
# snp_col=$(get_col "SNP")         # retrieve the index for column named "SNP"
# awk -v snp_col="$snp_col" ...

# --- get_file_delimiter() ---
# purpose: guess and return the file delimiter from the first line
# usage: get_file_delimiter "myfile.csv"
get_file_delimiter() {
    local file="$1"
    local header
    if [[ "$file" == *.gz ]]; then
        header=$(zcat "$file" | head -n 1)
    else
        header=$(head -n 1 "$file")
    fi

    # list of candidate delimiters (actual tab, not string '\t')
    local delimiters=(',' $'\t' ';' '|' ':' ' ')
    local max_count=0
    local guess=""
    local delim count

    for delim in "${delimiters[@]}"; do
        count=$(awk -F"$delim" '{print NF-1}' <<< "$header")
        if (( count > max_count )); then
            max_count=$count
            guess=$delim
        fi
    done

    # Only normalize spaces if the guess is a run of spaces and not a tab
    if [[ "$guess" =~ ^[[:space:]]+$ ]] && [[ "$guess" != $'\t' ]]; then
        guess=" "
    fi

    # Print the delimiter (tabs and spaces correctly)
    if [[ "$guess" == $'\t' ]]; then
        echo $'\t'
    elif [[ "$guess" == ' ' ]]; then
        echo " "
    else
        echo "$guess"
    fi
}

# --- parse_header() ---
# purpose: function to parse header and build associative array of column names (sets global col_indices)
# usage: parse_header "$file" "$delimiter"
parse_header() {
    local file="$1"
    local header
    if [[ "$file" == *.gz ]]; then
        header=$(zcat "$file" | head -n 1 || true)
    else
        header=$(head -n 1 "$file")
    fi

    local delimiter="$2"

    # make sure header provided
    if [[ -z "$header" ]]; then
        echo "Error: No header line provided to parse_header." >&2
        return 1
    fi

    # use delimiter
    IFS="$delimiter" read -r -a cols <<< "$header"

    # make an associative array between column names and indices, e.g. col_indices["snp"]=1, col_indices["chrom"]=2
    # sed line is just stripping spaces and other problematic characters in column name
    declare -gA col_indices
    for i in "${!cols[@]}"; do
        colname=${cols[$i]}
        colname=$(echo "$colname" | sed 's/^[`"'\'' ]*//;s/[`"'\'' ]*$//')
        col_indices[$colname]=$((i+1))
    done
    return 0
}

# --- get_col() ---
# purpose: takes column name-index associative array and prints out results and warnings
# usage: get_col "your_column_name"
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

# --- for debugging, section for standalone use ---
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
