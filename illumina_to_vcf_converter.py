"""
File: illumina_to_vcf_converter.py
Author: rezza.mozafari@gmail.com
Created: August 1, 2024

# Description
This script converts genotype data extracted from the Illumina format to VCF format,
considering the ILMN and Customer strand information.

# Usage
python illumina_to_vcf_converter.py --genotype_file <path_to_genotype_file> --marker_file <path_to_marker_file> --output_vcf <path_to_output_vcf> --chromosome <chromosome_number> --phased/--unphased

# Python Version
Developed and tested with Python 3.8.5

# Dependencies
- pandas==1.2.4
- numpy==1.20.1
- os
- time
- sys
- argparse
"""

import pandas as pd
import numpy as np
import os
import time
import sys
import argparse

def convert_strand(allele, ilmn_strand, customer_strand):
    """
    Convert the allele based on the strand information.
    
    Parameters:
    allele (str): The allele to convert (A, T, C, G).
    ilmn_strand (str): The strand information from the ILMN (Illumina) system.
    customer_strand (str): The strand information from the customer or reference.

    Returns:
    str: The converted allele based on the strand information.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    if ilmn_strand == customer_strand:
        return allele  # No conversion needed if strands match
    else:
        return complement.get(allele, '.')  # Convert allele if strands differ

def convert_genotype(genotype, ref, alt, phased):
    """
    Convert numeric genotype to VCF format.

    Parameters:
    genotype (array): The genotype data to convert.
    ref (str): The reference allele.
    alt (str): The alternate allele.
    phased (bool): Indicates if the genotype data is phased.

    Returns:
    array: The converted genotype in VCF format.
    """
    if phased:
        genotype_mapping = {
            '0': "{0}|{0}".format(alt if alt != '.' else '.'),
            '1': "{0}|{1}".format(ref if ref != '.' else '.', alt if alt != '.' else '.'),
            '2': "{0}|{0}".format(ref if ref != '.' else '.'),
            '5': './.',
            '-': './.'
        }
    else:
        genotype_mapping = {
            '0': "{0}/{0}".format(alt if alt != '.' else '.'),
            '1': "{0}/{1}".format(ref if ref != '.' else '.', alt if alt != '.' else '.'),
            '2': "{0}/{0}".format(ref if ref != '.' else '.'),
            '5': './.',
            '-': './.'
        }
    
    return np.vectorize(lambda x: genotype_mapping.get(x, './.'))(genotype)

def main(genotype_file, marker_file, output_vcf, chromosome, phased):
    try:
        print("Loading input files...")
        
        # Check and load the genotype file
        if os.path.exists(genotype_file):
            genotypes = pd.read_csv(genotype_file)
            print("Genotype file loaded successfully: %s" % genotype_file)
            print("Number of samples: %d" % len(genotypes))
        else:
            raise FileNotFoundError("Genotype file not found: %s" % genotype_file)

        # Check and load the marker file
        if os.path.exists(marker_file):
            map_file = pd.read_csv(marker_file, sep=";")
            print("Marker file loaded successfully: %s" % marker_file)
            print("Number of markers: %d" % len(map_file))
        else:
            raise FileNotFoundError("Marker file not found: %s" % marker_file)

        print("File loading completed successfully.")

        # Convert genotypes to an array for speed
        genotype_array = np.array([list(g) for g in genotypes['genotype']])

        print("Starting VCF file processing...")
        start_time = time.time()

        vcf_header = [
            "##fileformat=VCFv4.2",
            "##source=custom_script",
            "##reference=UMD3.1",
            '##FILTER=<ID=PASS,Description="All filters passed">',
            '##FILTER=<ID=LowQual,Description="Low quality SNP or missing allele information">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % "\t".join(genotypes['Matricola'])
        ]

        with open(output_vcf, "w") as vcf_file:
            # Write header
            for line in vcf_header:
                vcf_file.write(line + "\n")
            
            # Process markers
            total_markers = len(map_file)
            for index, marker in map_file.iterrows():
                if chromosome and marker['Chromosome'] != chromosome:
                    continue  # Skip markers not on the specified chromosome

                if index % 100 == 0:  # Update progress every 100 markers
                    elapsed_time = time.time() - start_time
                    estimated_total_time = (elapsed_time / (index + 1)) * total_markers
                    progress = f"\rProcessing marker {index} of {total_markers} ({index/total_markers*100:.2f}%). Estimated remaining time: {(estimated_total_time - elapsed_time) / 60:.2f} minutes"
                    sys.stdout.write(progress)
                    sys.stdout.flush()
                
                try:
                    chrom = marker['Chromosome']
                    pos = marker['Position']
                    snp_id = marker['Name']
                    snp = str(marker['SNP']) if pd.notnull(marker['SNP']) else ''
                    ilmn_strand = marker['ILMN_Strand']
                    customer_strand = marker['Customer_Strand']
                    
                    # Handle invalid SNPs using '.' instead of 'N'
                    if len(snp) != 2:
                        ref, alt = '.', '.'
                        filter_status = "LowQual"
                    else:
                        # Convert SNP alleles to customer strand
                        ref = convert_strand(snp[0], ilmn_strand, customer_strand)
                        alt = convert_strand(snp[1], ilmn_strand, customer_strand)
                        
                        filter_status = "PASS" if ref != '.' and alt != '.' else "LowQual"

                    genotype_values = convert_genotype(genotype_array[:, index], ref, alt, phased)

                    vcf_row = "{0}\t{1}\t{2}\t{3}\t{4}\t.\t{5}\t.\tGT\t{6}\n".format(
                        chrom, pos, snp_id, ref, alt, filter_status, "\t".join(genotype_values)
                    )
                    vcf_file.write(vcf_row)
                except KeyError as e:
                    print(f"\nError reading data for marker {index}: {e}")
                    continue
                except IndexError as e:
                    print(f"\nError accessing genotype data for marker {index}: {e}")
                    continue
                except Exception as e:
                    print(f"\nUnexpected error during processing of marker {index}: {e}")
                    continue

        print("\nVCF file created successfully: %s" % output_vcf)
        print("Total processing time: %.2f minutes" % ((time.time() - start_time) / 60))

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except pd.errors.EmptyDataError:
        print("Error: One of the input files is empty.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert Illumina genotype data to VCF format.')
    parser.add_argument('--genotype_file', type=str, required=True, help='Path to the genotype CSV file.')
    parser.add_argument('--marker_file', type=str, required=True, help='Path to the marker CSV file.')
    parser.add_argument('--output_vcf', type=str, required=True, help='Path to the output VCF file.')
    parser.add_argument('--chromosome', type=int, help='Chromosome number to filter data (optional).')
    parser.add_argument('--phased', action='store_true', help='Indicate if the genotype data is phased.')
    parser.add_argument('--unphased', action='store_true', help='Indicate if the genotype data is unphased.')

    args = parser.parse_args()
    
    # Ensure that either phased or unphased is specified
    if args.phased and args.unphased:
        raise ValueError("Please specify either --phased or --unphased, not both.")
    elif not args.phased and not args.unphased:
        raise ValueError("Please specify whether the input data is --phased or --unphased.")

    main(args.genotype_file, args.marker_file, args.output_vcf, args.chromosome, args.phased)