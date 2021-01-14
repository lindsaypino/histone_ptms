import sys
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import re
from Bio import SeqIO
from tqdm import tqdm

sys.stdout.write("Imported required packages successfully.\n")


# set the master table for modification mass shifts
MZSHIFT_DICT = {'[n-term PR]': '[+56]',
                '[AC]': '[+42]',
                '[PR]': '[+56]',
                '[ME1]': '[+14]',
                '[ME1+PR]': '[+70]',  # methyl+propionyl
                '[ME2]': '[+28]',
                '[ME3]': '[+42]',
                '[PH]': '[+80]',
                '[n-term ME3+PR]': '[+98.1]',
                '[n-term ME2+PR]': '[+84.1]',
                '[n-term ME1+PR]': '[+126.1]',
                '[n-term PR2]': '[+112.1]',
                '[n-term AC]': '[+98]',
                '[n-term AC+PR]': '[+84.1]',
                '[GGprop+PR]': '[+170.1]'}

# usage statement and input descriptions
parser = argparse.ArgumentParser(
    description='A post-processing script to decode Skyline\'s Peptide Modified Sequence value into a short-hand \
                    histone mark. Please note that the mass shifts (e.g. [+42] = ac, [+80] = ph) are hard-coded \
                    so if you have novel histone modifications in your Skyline document, you\'ll need to add them\
                    to this code!',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('fasta_file', type=str,
                    help='a file with protein names and protein sequences in fasta format (e.g. download from Uniprot')
parser.add_argument('skyline_groupcomp', type=str,
                    help='an export from Skyline\'s Group Comparison feature (View > Other Grids > Group Comparisons) \
                        which contains columns for Protein Name, Peptide Modified Sequence, Fold Change Result, and \
                        Adjusted P-Value. See user manual for tutorial on how to set this up in Skyline.')

# parse arguments from command line
args = parser.parse_args()
fasta_file = args.fasta_file
skyline_file = args.skyline_groupcomp

##
## read input files: FASTA and Skyline Export Report with Peptide Modified Sequences
##

# parse FASTA for protein name and sequence

protein_df = pd.DataFrame()  # Initialize a dataframe to store results
for protein in tqdm(SeqIO.parse(fasta_file, "fasta")):

    protein_sequence = str(protein.seq).upper()

    # cleave initial "start" methionine if present
    if protein_sequence[0] == "M":
        protein_sequence = protein_sequence[1:]

    # add this protein to the dataframe
    new_df = pd.DataFrame({'protein': protein.id,
                           'protein_sequence': protein_sequence}, index=[0])

    protein_df = protein_df.append(new_df)
protein_df = protein_df.drop_duplicates()

# read in Skyline Export Report with Peptide Modified Sequences
skyline_df = pd.read_csv(skyline_file)


##
## "decode" modified peptide sequences to biological histone marks
##

# remove propionylations ([+56], [+112.1]) which aren't biologically relevant here
skyline_df['new_pep_seq'] = skyline_df['Peptide Modified Sequence']
skyline_df['new_pep_seq'] = skyline_df['new_pep_seq'].str.replace(r'\[\+56\]', '')
skyline_df['new_pep_seq'] = skyline_df['new_pep_seq'].str.replace(r'\[\+112.1\]', '')

# decode each modified peptide sequence to its modified residue number and mod type
decode_df = pd.DataFrame()  # Initialize a dataframe to store results
for index, row in skyline_df.iterrows():
    peptide = row['Peptide']
    mod_seq = row['new_pep_seq']
    fc = row['Fold Change Result']
    pval = row['Adjusted P-Value']

    # find all protein matches for the peptide in case there are duplicate/non-unique peptides
    protein_match = list(protein_df[protein_df['protein_sequence'].str.contains(peptide)]['protein'])

    # map each peptide to its residue position in the protein and change mass shift to modification
    for protein in protein_match:
        sequence = protein_df[protein_df['protein'] == protein]['protein_sequence'][0]
        aa_index = sequence.find(peptide)  # get the amino acid index position

        # split the peptide sequence into residues including any [+nn] mod shift
        residue_list = re.sub(r"([A-Z])", r" \1", mod_seq).split()

        # build the [Residue][Index Position][modification] histone mark
        histone_mod = ''
        for i in range(len(residue_list)):
            if '[' in residue_list[i]:
                aa_pos = aa_index + i + 1  # have to +1 for indexing

                # match the mass shift to the mod and format to make it pretty
                mod = list(MZSHIFT_DICT.keys())[list(MZSHIFT_DICT.values()).index(residue_list[i][1:])]
                mod = re.sub('n-term ', '', mod)
                mod = re.sub(r'\+PR', '', mod)
                mod = mod.lower()

                # build [Residue][Index Position][modification] string
                new_mod = residue_list[i][0] + str(aa_pos) + mod
                histone_mod = histone_mod + new_mod

        # add this protein/peptide to the dataframe
        new_df = pd.DataFrame({'Protein Name': protein,
                               'Peptide Sequence': peptide,
                               'Peptide Modified Sequence': row['Peptide Modified Sequence'],
                               'histone mark': histone_mod,
                               'Fold Change Result': fc,
                               'Adjusted P-Value': pval}, index=[0])
        decode_df = decode_df.append(new_df)

decode_df = decode_df.drop_duplicates()

decode_df = decode_df.groupby(['Peptide Modified Sequence',
                               'histone mark',
                               'Fold Change Result',
                               'Adjusted P-Value'])['Protein Name'].apply(
    lambda x: ','.join(x)).reset_index()

decode_df.to_csv("D:/Penn/proj/collab_greer/data/greer_onlyhistlibrary_groupcomparison_mound-v-vegetative_decoded.csv",
                 index=False)
