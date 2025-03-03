# CRISPR Guide RNA Evaluation Tool

## Overview
This tool evaluates the suitability of guide RNAs for CRISPR knockout experiments in homo sapiens. The default example analyzes EGFR, but you can use it with any gene of interest.

## Installation

1. Download the zip file "Arjunsinh_Harer_Final_Project_RBIF_109"
2. Extract the contents
3. Save to your desktop
4. Navigate to the folder in your terminal

## Basic Usage

To see an example analysis (EGFR knockout), run:
```bash
python off_target.py
```

## Analyzing Your Own Gene

### Step 1: Download Gene Data
Run the following command:
```bash
python download_fasta.py "Gene Name" "Your Email"
```

Example:
```bash
python download_fasta.py EGFR arjunsinhharer@brandeis.edu
```

### Step 2: Get Guide RNA Sequences

1. Visit the [Synthego Knockout Guide Design Tool](https://design.synthego.com/)
2. Set genome type to "Homo sapiens"
3. Enter your gene of interest
4. Click "Search"
5. Copy the recommended guide sequences

### Step 3: Prepare Guide RNA File

1. Remove any existing .txt files from the project directory
2. Create a new text file named `guides_{GeneName}.txt`
   - Example: `guides_EGFR.txt`
3. Paste the guide sequences from Synthego
   - One sequence per line
4. Save the file in the project directory

### Step 4: Run Analysis
```bash
python off_target.py
```

## Important Notes

- The script requires exactly one .txt file in the directory to run properly
- File naming convention must be followed: `guides_{GeneName}.txt`
- The FASTA file for your gene must be downloaded before analysis

## Support

If you encounter any issues, please reach out

- 📧 Email: arjunsinhharer@brandeis.edu

## Author

Arjunsinh Harer  
Brandeis University  
RBIF 109 Final Project
