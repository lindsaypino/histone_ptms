import os
import sys
import argparse
import pandas as pd

# usage statement and input descriptions
parser = argparse.ArgumentParser(
    description='A post-processing script to merge files from Skyline\'s Group Comparison report format.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('annotation_file', type=str,
                    help='a file with each of the filenames and the annotation')
parser.add_argument('--output_path', default=os.getcwd(), type=str,
                    help='specify an output path for the merged result file')


# parse arguments from command line
args = parser.parse_args()
annotation_file = args.annotation_file
output_path = args.output_path

cwd = os.getcwd()

# iterate over the filenames, add a sample group column, and merge them together
annotation_df = pd.read_csv(annotation_file)
files2merge = annotation_df['filename']
merge_df = pd.DataFrame()
for file in files2merge:
    temp_df = pd.read_csv(os.path.join(cwd, file))
    annotation = annotation_df[annotation_df['filename'] == file]['annotation'].iloc[0]
    temp_df['samplegroup'] = [annotation] * len(temp_df)
    merge_df = merge_df.append(temp_df)

merge_df.to_csv(path_or_buf=os.path.join(output_path, 'merged_skyline_groupcomparisons.csv'),
                index=False)

sys.stdout.write("Finished merging Skyline Group Comparison exports.\n")