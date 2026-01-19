#!/usr/bin/env bash
set -euo pipefail

# User-configurable inputs
vcf_file="/Users/k.zoltowska/Documents/General_scripts/NEUROX_filtered.vcf"
output_prefix="$(basename "$vcf_file" .vcf)"

gene_list="LRRK2,PINK1,PARKN,VPS35,PARK7"
rsid_csv="/Users/k.zoltowska/Documents/General_scripts/PD_variants.csv"   # CSV must contain column 'rsid'

subset_for_region=false
region_chr="12"
region_from=40657638
region_to=40713560

genome_build="hg19"   # ANNOVAR uses hg19/hg38 notation

# Paths for ANNOVAR reference database
annovar_db="/Users/k.zoltowska/Documents/software/annovar/humandb"   # update to your humandb path

vcf_tmp="vcf_tmp.vcf"
vcf_clean="vcf_clean.vcf.gz"

# Optional region subset

if [ "$subset_for_region" = true ]; then
    vcftools \
      --vcf "$vcf_file" \
      --chr "$region_chr" \
      --from-bp "$region_from" \
      --to-bp "$region_to" \
      --recode \
      --recode-INFO-all \
      --out "${vcf_tmp%.vcf}"
    vcf_input="${vcf_tmp%.vcf}.recode.vcf"
else
    vcf_input="$vcf_file"
fi

# Annotate with ANNOVAR
annotated_prefix="${output_prefix}_annovar"

echo "Annotating with ANNOVAR..."
/Users/k.zoltowska/Documents/software/annovar/table_annovar.pl "$vcf_input" "$annovar_db" \
    -buildver "$genome_build" \
    -out "$annotated_prefix" \
    -remove \
    -protocol refGene \
    -operation g \
    -nastring "." \
    -vcfinput

annotated_vcf="${annotated_prefix}.hg19_multianno.vcf"
# Compress + index
bgzip -f "$annotated_vcf"
bcftools index "${annotated_vcf}.gz"

# Filter by gene list

echo "Filtering annotated VCF by gene list..."
gene_vcf="${output_prefix}_gene_filtered.vcf"

gene_expr=$(echo "$gene_list" | awk -F',' '{
  for(i=1;i<=NF;i++){
    printf "\''%s'\''", $i
    if(i<NF){printf "|"}
  }
}')

echo ${gene_expr}

# Using bcftools + ANNOVAR INFO field
bcftools view \
  -i "INFO/Gene.refGene ~ \"(^|,)(${gene_expr})(,|$)\"" \
  "${annotated_vcf}.gz" -Oz -o "${gene_vcf}.gz"
bcftools index "${gene_vcf}.gz"

# Filter by rsID list using bcftools
echo "Filtering annotated VCF by rsID list..."
rsid_vcf="${output_prefix}_rsid_filtered.vcf"

tail -n +2 "$rsid_csv" | cut -d',' -f2 | sort -u > rsids.txt

# Extract variants whose IDs contain any rsID substring
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "${annotated_vcf}.gz" | \
# this allows partial matches (good as we include Exm_rs12121 but bad as we include rs123 for rs12)
awk 'BEGIN{
    while(getline < "rsids.txt") rsid[$1]=1
}
{
    for(id in rsid)
        if(index($3,id))
            print $1"\t"$2"\t"$3"\t"$4"\t"$5 
}' > rsid_positions.txt

# Convert to proper region file for bcftools
awk '{start=$2-1; if(start<0) start=0; end=$2; print $1"\t"start"\t"end}' rsid_positions.txt > rsid_regions.txt

# Extract variants from annotated VCF
bcftools view -R rsid_regions.txt "${annotated_vcf}.gz" -Oz -o "${rsid_vcf}.gz"
bcftools index "${rsid_vcf}.gz"

# Cleanup
rm rsids.txt rsid_positions.txt rsid_regions.txt

# Combine gene + rsID filtered VCFs (union, no duplicates)
combined_vcf="${output_prefix}_protein_coding_combined.vcf.gz"

# D drops duplicates, a allows same positions
bcftools concat -a -D "${gene_vcf}.gz" "${rsid_vcf}.gz" | bcftools sort -Oz -o "$combined_vcf"
bcftools index "$combined_vcf"

# Generate genotype matrix
genotype_prefix="${output_prefix}_genotype_matrix"
vcftools --gzvcf "$combined_vcf" --012 --out "$genotype_prefix"

# Get the labels for the columns in the .tsv file
bcftools query -f '%ID\t%CHROM\t%POS\t%INFO\n' "$combined_vcf" | \
awk -F'\t' '
{
    # Use rsID if present, otherwise CHROM_POS
    if ($1 != ".") {
        vid = $1
    } else {
        vid = $2 "_" $3
    }

    # Replace characters that break headers (optional but recommended)
    gsub(/[[:space:]]+/, "", $4)

    # Print: id_chrom_pos_annotations
    print vid "_" $2 "_" $3 "_" $4
}' > "${genotype_prefix}_variant_labels.txt"


# Create TSV
echo "Creating tsv"
awk '{for(i=2;i<=NF;i++){if($i==-1)$i="NA"; printf "%s%s",$i,(i==NF?ORS:OFS)}}' \
OFS='\t' "${genotype_prefix}.012" > "${genotype_prefix}.012.clean"

echo -e "Sample\t$(paste -sd'\t' "${genotype_prefix}_variant_labels.txt")" > "${genotype_prefix}.tsv"
paste -d'\t' "${genotype_prefix}.012.indv" "${genotype_prefix}.012.clean" >> "${genotype_prefix}.tsv"


# Cleanup

echo "All done!"
