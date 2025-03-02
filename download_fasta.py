# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 17:33:52 2024

@author: Arjunsinh Harer
"""

from Bio import Entrez
import sys

def download_fasta(gene_name, email, database="nucleotide", max_results=1):
    """
    Download FASTA files for a given gene from NCBI, limited to Homo sapiens.

    Parameters:
    - gene_name: The name of the gene to search for.
    - email: Your email address (required by NCBI).
    - database: The NCBI database to search in. Default is 'nucleotide'.
    - max_results: Maximum number of search results to fetch. Default is 1.
    """
    Entrez.email = email  # Always tell NCBI who you are
    # Append '[ORGN]' to the gene name to specify organism
    search_handle = Entrez.esearch(db=database, term=f"{gene_name}[Gene Name] AND Homo sapiens[ORGN]", retmax=max_results)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    ids = search_results['IdList']
    if not ids:
        print(f"No results found for gene: {gene_name} in Homo sapiens")
        return

    for i, gene_id in enumerate(ids, start=1):
        fetch_handle = Entrez.efetch(db=database, id=gene_id, rettype="fasta", retmode="text")
        fasta_data = fetch_handle.read()
        fetch_handle.close()

        # Save the FASTA data to a file
        filename = f"{gene_name}_homo_sapiens_{i}.fasta"
        with open(filename, 'w') as file:
            file.write(fasta_data)
        print(f"Downloaded {filename}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <gene_name> <email>")
    else:
        gene_name = sys.argv[1]
        email = sys.argv[2]
        download_fasta(gene_name, email)
