This repository contains a Python script designed for bioinformatics tasks, specifically to extract and format features from GFF3 files. 
The goal is to provide an easy-to-use tool for parsing genomic data and retrieving sequences in FASTA format.

Features

	•	Extracts specific features from a GFF3 file based on their type and attribute values.
	•	Handles FASTA sections in GFF3 files, enabling sequence extraction.
	•	Supports reverse complement generation for features on the negative strand.
	•	Formats output sequences in FASTA format with customizable line lengths.

Requirements

	•	Python 3.6 or higher.

Installation

Clone the repository to your local machine:

	•	git clone https://github.com/cgregg1216/Bioinformatic_Portfolio.git
	•	cd Bioinformatic_Portfolio

Usage

Run the script with the following arguments:

	•	python3 export_gff3_feature.py --source_gff <path_to_gff3_file> --type <feature_type> --attribute <attribute_key> --value <attribute_value>

Example

To extract a gene with the ID YAL069W from a GFF3 file named example.gff3:

	•	python3 export_gff3_feature.py --source_gff example.gff3 --type gene --attribute ID --value YAL069W

The output will be displayed in FASTA format:

>gene:ID:YAL069W
ATGCTAGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTAGCTAGCTAGC

Arguments

	•	--source_gff: Path to the GFF3 file.
	•	--type: Type of feature to extract (e.g., gene, mRNA).
	•	--attribute: Attribute key to match (e.g., ID, Name).
	•	--value: Attribute value to match (e.g., YAL069W).

Functions

	•	parse_gff3: Parses the GFF3 file and retrieves features and sequence data.
	•	get_sequence: Extracts the sequence for a feature using its coordinates and strand.
	•	reverse_complement: Generates the reverse complement of a DNA sequence.
	•	format_fasta: Formats sequences into FASTA format.

Notes

	•	The script alerts you if no matching features are found or if multiple matches exist.
	•	Proper error handling is implemented for invalid inputs or incomplete data.

Contribution

Feel free to fork the repository and submit pull requests for improvements or bug fixes.

License

This project is licensed under the MIT License. See the LICENSE file for details.
