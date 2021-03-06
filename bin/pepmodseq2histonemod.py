import sys
import os
import argparse
import pandas as pd
import re
from Bio import SeqIO


# set the master table for modification mass shifts
MZSHIFT_DICT = {'[n-term PR]': '[+56]',
                '[AC]': '[+42]',
                '[PR]': '[+56]',
                '[ME1]': '[+14]',
                '[ME1+PR]': '[+70]',  # methyl+propionyl
                '[ME2]': '[+28]',
                '[ME3]': '[+42]',
                '[PH]': '[+80]',
                '[PH+PR]': '[+136]',
                '[n-term ME3+PR]': '[+98.1]',
                '[n-term ME2+PR]': '[+84.1]',
                '[n-term ME1+PR]': '[+126.1]',
                '[n-term PR2]': '[+112.1]',
                '[n-term AC]': '[+98]',
                '[n-term AC+PR]': '[+84.1]',
                '[GG+PR]': '[+170.1]',
                '[GG]': '[+114]',
                '[OX]': '[+16]',
                '[D3-AC]': '[+3]',
                '[ME2+D3ac]': '[+9.1]'}

# usage statement and input descriptions
parser = argparse.ArgumentParser(
    description='A post-processing script to decode Skyline\'s Peptide Modified Sequence value into a short-hand \
                    histone mark. Please note that the mass shifts (e.g. [+42] = ac, [+80] = ph) are hard-coded \
                    so if you have novel histone modifications in your Skyline document, you\'ll need to add them\
                    to this code! This code is intended for use with Skyline\'s Document Grid>Peptides report format.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('fasta_file', type=str,
                    help='a file with protein names and protein sequences in fasta format. The protein names must \
                         match the Protein in the Skyline document (use the same fasta!)')
parser.add_argument('docgrid_peptides', type=str,
                    help='an export from Skyline\'s "Peptides" document grid report (View > Document Grid > Peptides) \
                        which contains the minimal required columns for Protein Name and Peptide Modified Sequence. \
                        See user manual for tutorial for more details on how to set this up in Skyline.')
parser.add_argument('--output_path', default=os.getcwd(), type=str,
                    help='specify an output path for the decoded result file')

# parse arguments from command line
args = parser.parse_args()
fasta_file = args.fasta_file
skyline_file = args.docgrid_peptides
output_dir = args.output_path


##
## read input files: FASTA and Skyline Export Report with Peptide Modified Sequences
##

# parse FASTA for protein name and sequence

protein_df = pd.DataFrame()  # Initialize a dataframe to store results
for protein in SeqIO.parse(fasta_file, "fasta"):

    protein_sequence = str(protein.seq).upper()

    # add some protein description parsing for non-model organisms that don't have useful gene names
    if "|" and "OS=" in str(protein.description):
        protein_desc = re.sub(r'([a-z]+\|[A-Z][0-9A-Z]+\.?\d*\|[0-9A-Z]+_[A-Z]+ )', '', str(protein.description)).split(" OS=", 1)[0]

    # cleave initial "start" methionine if present
    if protein_sequence[0] == "M":
        protein_sequence = protein_sequence[1:]

    # add this protein to the dataframe
    new_df = pd.DataFrame({'protein': protein.id,
                           'protein_sequence': protein_sequence,
                           'protein_desc': protein_desc}, index=[0])

    protein_df = protein_df.append(new_df)
protein_df = protein_df.drop_duplicates()

# read in Skyline Export Report with Peptide Modified Sequences
skyline_df = pd.read_csv(skyline_file)

# check that the Skyline Export Report format has the minimum requirements
if 'Peptide' in skyline_df.columns:
    pass
elif 'Peptide Sequence' in skyline_df.columns:
    skyline_df = skyline_df.rename(columns={"Peptide Sequence": "Peptide"})
else:
    sys.exit('ERROR: Skyline export file must include Peptide Sequence column.\n')
if 'Peptide Modified Sequence' in skyline_df.columns:
    pass
else:
    sys.exit('ERROR: Skyline export file must include Peptide Modified Sequence column.\n')
if 'Protein' in skyline_df.columns:
    pass
elif 'Protein Name' in skyline_df.columns:
    skyline_df = skyline_df.rename(columns={"Protein Name": "Protein"})
else:
    sys.exit('ERROR: Skyline export file must include Protein column.\n')

sys.stdout.write("Successfully imported data, decoding modified peptide sequences now.\n")


##
## "decode" modified peptide sequences to biological histone marks
##

skyline_df['new_pep_seq'] = skyline_df['Peptide Modified Sequence']

# decode each modified peptide sequence to its modified residue number and mod type
decode_df = pd.DataFrame()  # Initialize a dataframe to store results
for index, row in skyline_df.iterrows():

    peptide = row['Peptide']
    mod_seq = row['new_pep_seq']

    # find all protein matches for the peptide in case there are duplicate/non-unique peptides
    protein_matches = []
    protein_matches = list(protein_df[protein_df['protein_sequence'].str.contains(peptide)]['protein'])

    # map each peptide to its residue position in the protein and change mass shift to modification
    for protein in protein_matches:
        sequence = protein_df[protein_df['protein'] == protein]['protein_sequence'][0]
        aa_index = sequence.find(peptide)  # get the amino acid index position

        # split the peptide sequence into residues including any [+nn] mod shift
        residue_list = re.sub(r"([A-Z])", r" \1", mod_seq).split()

        # build the [Residue][Index Position][modification] histone mark
        histone_mod = ''
        for i in range(len(residue_list)):
            if '[' in residue_list[i] and '56' not in residue_list[i]:
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
        new_df = row
        new_df['histone mark'] = histone_mod
        new_df['protein_match'] = protein
        new_df['protein_desc'] = protein_df[protein_df['protein'] == protein]['protein_desc'][0]
        #print(new_df)

        decode_df = decode_df.append(new_df)

decode_df = decode_df[['new_pep_seq', 'protein_match', 'protein_desc', 'histone mark']]
decode_df.drop_duplicates(inplace=True)

allbut = list(decode_df.columns)
allbut.remove('protein_match')
allbut.remove('protein_desc')

# make protein "groupings" for each non-unique histone mark
concat_names = decode_df.groupby(allbut)['protein_desc'].apply(lambda x: ','.join(x)).reset_index(); print(concat_names.head())
decode_df = decode_df.groupby(allbut)['protein_match'].apply(lambda x: ','.join(x)).reset_index(); #print(decode_df.head())
decode_df['protein_desc'] = concat_names['protein_desc']
decode_df.head()

#sys.exit()

#decode_df.drop(['Protein'], axis=1, inplace=True)
decode_df.rename(columns = {'new_pep_seq': 'Peptide Modified Sequence',
                            'protein_match': 'Protein Group',
                            'protein_desc': 'Protein Description'},
                 inplace = True)

decode_df = decode_df[['Peptide Modified Sequence', 'Protein Group', 'Protein Description', 'histone mark']]

decode_df.drop_duplicates(inplace=True)

# TODO fix this parsing to get just the filename
out_file = os.path.splitext(skyline_file)[0]; print(out_file)

decode_df.to_csv(path_or_buf=os.path.join(output_dir, (out_file+'_decoded.csv')), index=False)

sys.stdout.write("Finished decoding modified peptide sequences.\n")
