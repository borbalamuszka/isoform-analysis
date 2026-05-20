import pandas as pd
import sys
import os

def convert_matrix(input_path, sep_out='\t'):
    # Generate output path by replacing .txt with .tsv
    output_path = input_path.rsplit('.txt', 1)[0] + '.tsv'
    
    # Read the input file (assuming tab-separated)
    df = pd.read_csv(input_path, delimiter='\t', index_col=0)
    # Write to output file
    df.to_csv(output_path, sep=sep_out)
    print(f"Converted {input_path} to {output_path}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python convert_to_tsv.py <input_file.txt>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    if not os.path.exists(input_file):
        print(f"Error: File {input_file} does not exist")
        sys.exit(1)
    
    if not input_file.endswith('.txt'):
        print("Error: Input file must be a .txt file")
        sys.exit(1)
    
    convert_matrix(input_file)