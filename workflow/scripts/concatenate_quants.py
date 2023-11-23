
import argparse

def open_quant_files(input_file, delimiter="\t"):
    with open(input_file, "r") as f:
        header = next(f, None)

        if header is None:
            return {}

        counts = {}
        lines = f.readlines()

        for line in lines:
            key, values = line.strip().split(delimiter, 1)
            counts[key] = values

    return counts

def concatenate_files(input_list):
    # List to store dictionaries from each file
    dictionaries = []

    # Open and read each file, store the dictionary in the list
    for file_path in open(input_list,"r"):
        counts = open_quant_files(file_path.strip('\n'))
        dictionaries.append(counts)

    # Compare keys of dictionaries
    if len(dictionaries) > 1:
        
        
        common_keys = set(dictionaries[0].keys())

        for counts in dictionaries[1:]:
            common_keys &= set(counts.keys())

        # Combine values for common keys
        combined_values = {}
        for key in common_keys:
            values_list = [counts[key] for counts in dictionaries]
            combined_values[key] = '\t'.join(values_list)

        print("Combined values for common keys:")
        for key, value in combined_values.items():
            print(f"{key}: {value}")
    else:
        print("Not enough files for comparison.")

# Example usage


## Implementation ##
def parse_arguments():
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating gff3 files.')
    parser.add_argument('-q','--quants', required=True, type=str, help='List of quants files you want to concatenate')
    parser.add_argument('-o','--out', required=False , type=str, help='Output file')
    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        print("Error: No arguments provided.")

    return parser.parse_args()


def main():
    
    parser = argparse.ArgumentParser(description='Concatenate quant files')
    args = parse_arguments()

    print ( concatenate_files(args.quants) )

    
if __name__ == "__main__":
    main()
    