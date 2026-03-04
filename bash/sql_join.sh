#!/bin/bash

# Inputs: 
#    trait_files - an array of TSV files, each with columns chrom, snp, REF, ALT, plus other trait columns
# Output:
#    temp.tsv - fully outer-joined result

DB="my_gwas.db"
FINAL_TABLE="merged_final"

# Remove previous DB if exists
rm -f "$DB"

# 1. Import all files into SQLite tables
#echo "Importing files into SQLite..."

#index=1
#declare -a tables

#for file in "${trait_files[@]}"; do
#    # Create a unique table name per file (no dashes or dots)
#    table="trait_${index}"
#    tables+=("$table")
#    sqlite3 "$DB" <<EOF
#DROP TABLE IF EXISTS $table;
#.mode tabs
#.import $file $table
#EOF
#    ((index+=1))
#done

#echo "${#tables[@]} tables imported."

# 2. Sequential FULL OUTER JOIN using SQL
echo "Joining tables..."

# Start with first table
sqlite3 "$DB" <<EOF
DROP TABLE IF EXISTS merged_1;
CREATE TABLE merged_1 AS SELECT * FROM ${tables[0]};
EOF

for ((i=1; i<${#tables[@]}; i++)); do
    prev="merged_$i"
    next="merged_$((i+1))"
    curr="${tables[$i]}"

    sqlite3 "$DB" <<EOF
DROP TABLE IF EXISTS ${next};

CREATE TABLE ${next} AS
SELECT 
    COALESCE(p.chrom, c.chrom) AS chrom,
    COALESCE(p.snp,   c.snp)   AS snp,
    COALESCE(p.REF,   c.REF)   AS REF,
    COALESCE(p.ALT,   c.ALT)   AS ALT,
    p.*, c.*
FROM ${prev} p
LEFT JOIN ${curr} c
    ON  p.chrom = c.chrom
    AND p.snp   = c.snp
    AND p.REF   = c.REF
    AND p.ALT   = c.ALT

UNION

SELECT 
    COALESCE(p.chrom, c.chrom) AS chrom,
    COALESCE(p.snp,   c.snp)   AS snp,
    COALESCE(p.REF,   c.REF)   AS REF,
    COALESCE(p.ALT,   c.ALT)   AS ALT,
    p.*, c.*
FROM ${curr} c
LEFT JOIN ${prev} p
    ON  p.chrom = c.chrom
    AND p.snp   = c.snp
    AND p.REF   = c.REF
    AND p.ALT   = c.ALT
WHERE p.chrom IS NULL;
EOF
done

# The final merged table will be named merged_N
final="merged_${#tables[@]}"

# 3. Export merged result to temp.tsv
echo "Exporting final result to temp.tsv..."

sqlite3 "$DB" <<EOF
.headers on
.mode tabs
SELECT * FROM $final;
EOF > temp.tsv

echo "Done! Output in temp.tsv"
