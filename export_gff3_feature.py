#!/usr/bin/env python3

import argparse

def parse_gff3(gff3_file, feature_type, attribute_key, attribute_value):
    """
    Parse the GFF3 file and return the coordinates and strand of the specified feature.

    Args:
    gff3_file (str): Path to the GFF3 file.
    feature_type (str): Type of feature to look for (e.g., gene).
    attribute_key (str): Key of the attribute to match (e.g., ID).
    attribute_value (str): Value of the attribute to match (e.g., YAL069W).

    Returns:
    tuple: A tuple containing the features list and sequence data dictionary.
    """
    features = []  # List to store features that match the criteria
    fasta_start = False  # Flag to indicate the start of the FASTA section
    sequence_data = {}  # Dictionary to store sequence data
    current_sequence_id = None  # Variable to store the current sequence ID

    with open(gff3_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('##FASTA'):
                # Indicate that the FASTA section has started
                fasta_start = True
                continue
            if fasta_start:
                # Collect FASTA sequences
                if line.startswith('>'):
                    # Start of a new sequence
                    current_sequence_id = line[1:].split()[0]
                    sequence_data[current_sequence_id] = []
                else:
                    # Add sequence data to the current sequence ID
                    sequence_data[current_sequence_id].append(line)
            else:
                if line.startswith('#'):
                    # Skip comment lines
                    continue
                parts = line.split('\t')
                if len(parts) != 9:
                    # Skip lines that do not have 9 columns
                    continue
                seqid, source, ftype, start, end, score, strand, phase, attributes = parts
                if ftype != feature_type:
                    # Skip lines that do not match the feature type
                    continue
                # Parse the attributes column into a dictionary
                attr_dict = {key: value for key, value in [attr.split('=') for attr in attributes.split(';')]}
                if attr_dict.get(attribute_key) == attribute_value:
                    # If the attribute matches the specified key and value, add to features list
                    features.append((seqid, int(start), int(end), strand))

    return features, sequence_data

def get_sequence(sequence_data, seqid, start, end, strand):
    """
    Extract the sequence from the FASTA data based on the coordinates and strand.

    Args:
    sequence_data (dict): Dictionary containing sequence data.
    seqid (str): Sequence ID.
    start (int): Start coordinate of the feature.
    end (int): End coordinate of the feature.
    strand (str): Strand of the feature.

    Returns:
    str: The extracted sequence.
    """
    # Concatenate the sequence lines and extract the relevant portion
    sequence = ''.join(sequence_data[seqid])[start-1:end]
    if strand == '-':
        # If the strand is negative, return the reverse complement
        sequence = reverse_complement(sequence)
    return sequence

def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.

    Args:
    seq (str): DNA sequence.

    Returns:
    str: Reverse complement of the DNA sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def format_fasta(header, sequence, line_length=60):
    """
    Format the sequence into FASTA format with the specified line length.

    Args:
    header (str): Header for the FASTA format.
    sequence (str): DNA sequence.
    line_length (int): Number of characters per line in the FASTA format.

    Returns:
    str: Formatted FASTA string.
    """
    return f">{header}\n" + '\n'.join(sequence[i:i+line_length] for i in range(0, len(sequence), line_length))

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Export specific features from a GFF3 file.")
    parser.add_argument('--source_gff', required=True, help="Path to the GFF3 file")
    parser.add_argument('--type', required=True, help="Feature type (e.g., gene, mRNA)")
    parser.add_argument('--attribute', required=True, help="Attribute key (e.g., ID, Name)")
    parser.add_argument('--value', required=True, help="Attribute value (e.g., YAR003W)")

    # Parse arguments
    args = parser.parse_args()

    # Parse the GFF3 file for the specified feature
    features, sequence_data = parse_gff3(args.source_gff, args.type, args.attribute, args.value)
    if not features:
        # No matching features found
        print(f"No features found for {args.type}:{args.attribute}={args.value}")
        return
    if len(features) > 1:
        # More than one matching feature found
        print(f"Warning: More than one feature found for {args.type}:{args.attribute}={args.value}")
        for feature in features:
            seqid, start, end, strand = feature
            sequence = get_sequence(sequence_data, seqid, start, end, strand)
            header = f"{args.type}:{args.attribute}:{args.value}"
            print(format_fasta(header, sequence))
    else:
        # Exactly one matching feature found
        seqid, start, end, strand = features[0]
        sequence = get_sequence(sequence_data, seqid, start, end, strand)
        header = f"{args.type}:{args.attribute}:{args.value}"
        print(format_fasta(header, sequence))

if __name__ == "__main__":
    main()
