# Illumina Final Report txt file to VCF Converter

## Description
This script converts genotype data extracted from the Illumina Final Report txt file format to VCF format, considering the ILMN and Customer strand information. It is designed for use in genomic data analysis and is compatible with bioinformatics tools like PLINK and VCFtools.

## Features
- Converts Illumina Final Report file format to VCF format.
- Supports chromosome targeting for data selection.
- Handles strand information for accurate allele representation.
- Supports both phased and unphased genotype data.

## Sample Input and Output

### Sample Input Files

**Genotype Data**
```csv
Matricola,genotype
Sample1,0112010112502102
Sample2,0022010212002102
Sample3,5012010102002212
```

**Marker Data**
```csv
ID;Name;Chromosome;Position;GenTrain_Score;SNP;ILMN_Strand;Customer_Strand;NormID
1;ARS-BFGL-NGS-100024;19;28759876;0.7377;TC;B;T;0
2;ARS-BFGL-NGS-100098;19;29316818;0.7745;AG;T;B;0
3;ARS-BFGL-NGS-100175;19;46661667;0.782;AG;T;B;0
```

### Sample Output File

**Expected VCF Output**
```vcf
##fileformat=VCFv4.2
##source=custom_script
##reference=UMD3.1
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality SNP or missing allele information">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	01IT021001544317	01IT022990038158
19	29128377	ARS-BFGL-NGS-100024	A	G	.	PASS	.	GT	A/G	G/G	G/G
19	29688117	ARS-BFGL-NGS-100098	T	C	.	PASS	.	GT	C/C	T/C	C/C
19	47318036	ARS-BFGL-NGS-100175	T	C	.	PASS	.	GT	T/C	T/T	T/T
```

## Usage
To run the converter with the sample files, use the following command:

## Handling Phased and Unphased Data
The script is designed to handle both phased and unphased genotype data, but it requires the user to specify the format of the input data:

- **Input Format**: The input genotype data in `genotype.csv` does not explicitly indicate whether the data is phased or unphased. Therefore, it is essential for the user to know the format of their data before running the script.

- **User Option**: When running the script, the user must specify whether the input data is phased or unphased using the `--phased` argument:
  - Use `--phased` if the genotype data is in phased format (e.g., `0|1`).
  - Use `--unphased` if the genotype data is in unphased format (e.g., `0/1`).

The output VCF will reflect the genotype format used in the input data, ensuring compatibility with downstream analysis tools.

### Example Command
To run the converter, specifying the input format, use the following command:

```bash
python PATH/TO/illumina_to_vcf_converter.py --genotype_file PATH/TO/genotype.csv --marker_file PATH/TO/marker.csv --output_vcf PATH/TO/output.vcf --phased
```

## Requirements
Make sure to install the required dependencies listed in `requirements.txt` before running the script.

## Author
Reza Mozafari (rezza.mozafari@gmail.com)