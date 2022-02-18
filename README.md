# genofunc

A compiled tool for HIV sequence processing including referencing, annotation, feature extraction, filtering and merging.


### Fasta and/or Metadata Manipulation

Mini tools for manipulation of fasta and/or metadata files. Contains functions such as merging multiple fasta files and metadata files, filtering duplicates, encoding fasta file sequence IDs etc.


#### concatenate_fasta

Description: Concatenate all fasta file with the same gene regions into one sequence fasta file. 

Options:
    --in-prefix: Input prefix within a directory for specific fasta files to be read (Required)
    --gene: Gene regions to concatenate by which fasta file names are included in (Required)
    --out-dir: Output directory after concatenating all fasta files with the same gene name (Default: ./)
    --log-file: Output additional information from the program (default: stdout)

e.g. genofunc concatenate_fasta --in-prefix test/data/sequences/gene_fasta/ --gene pol gag --out-dir test/data/output/

#### extract_metadata

Description: Extract relevent metadata based on index fields from metadata file using fasta files. If field is empty, it will not be extracted and flagged in log file.

Options:
    --in-fasta: Input fasta file (Required)
    --in-metadata: Input metadata file (Required)
    --column: column name for extraction (Required)
    --rename: column name for re-naming (Default: None)
    --id-column: column name containing sequence ID (Required)
    --out-fasta: Output fasta file (Default: filtered.fasta)
    --out-metadata: Output metadata file (Default: extracted_metadata.csv)
    --log-file: Output additional information from the program (default: stdout)

e.g. genofunc extract_metadata --in-fasta test/data/sequences/seqB.fasta --in-metadata test/data/metadata/metadataB.csv --column country --id-column strain --out-fasta test/data/result/extract.fasta --out-metadata test/data/output/extract_metadata.csv

#### filter_fasta

Description: Filter fasta file based on minimum sequence length.

Options:
    --in-dir: Input directory (Required)
    --in-metadata: Input metadata for filtering alongside sequence file (Required)
    --genes: Which genes fasta files are needed for filtering (Required)
    --min-length: Minimum base coverage (Non-gaps) for sequences to be filtered if under this threshold (Required)
    --symmetric: Require all gene regions to be available for the same sequence (Default: False)
    --out-dir: Output directory after filtering (Default: ./)

e.g. genofunc filter_fasta --in-dir test/data/sequences/gene_fasta/ --in-metadata test/data/metadata/metadataB.csv --genes pol gag --min-length 0.95 --symmetric --out-dir test/data/output/

#### gene_concatenator

Description: Concatenate fasta file based on similar sequence names of multiple genomic regions. 

Options:
    --in-fasta: Multiple fasta files for concatenation (Required)
    --out-fasta: Output fasta file (Default: concatenate.fasta)

e.g. genofunc gene_concatenator --in-fasta test/data/sequences/gene_fasta/pol.fasta test/data/sequences/gene_fasta/gag.fasta --out-fasta test/data/output/gag_pol.fasta

#### merge

Description: Merges two or more fasta files avoiding duplicates based on matches to metadata

Options:
    --in-fasta: List of fasta files with spaces in between. At least two fasta files must be inserted here. (Required)
    --in-metadata: list of matching metadata file with same naming convention as fasta file. (Required)
    --index-field: The column ID with matching sequence IDs with fasta file (Required)
    --out-dir: Output merged sequence and metadata file from multiple inputs (Default: ./).

e.g. genofunc merge --in-fasta test/data/sequences/seqA.fasta test/data/sequences/seqB.fasta --in-metadata test/data/metadata/metadataA.csv test/data/metadata/metadataB.csv --index-field strain --out-dir test/data/output/

#### name_splitter

Description: Split the sequence name into metadata based on piping character.

Options:
    --in-fasta: Input fasta file (Required)
    --pipe: Input character the fasta sequence name is split based on (Required)
    --header: Header for the output metadata table (Default: "")
    --out-metadata: Output metadata file (Default: metadata.csv)

e.g. genofunc name_splitter --in-fasta test/data/sequences/seq_pipe.fasta --pipe | --out-metadata test/data/output/output_metadata.csv

#### rename_fasta

Description: Renaming fasta sequence names based on character splits.

Options:
    --in-fasta: Input fasta file (Required)
    --pipe: Input character the fasta sequence name is split based on (Required)
    --out-fasta: Output fasta file (Default: cleaned_sequences.fasta)

e.g. genofunc rename_fasta --in-fasta test/data/sequences/seq_pipe.fasta --pipe | --out-fasta test/data/output/depipe.fasta

#### strain_encoder

Description: Encoded strain id into non-defining ids. 

Options:
    --in-fasta: Input folder containing sequence files in fasta format (Required)
    --in-metadata: Input metadata corresponding to sequence files (Required)
    --encoding-column: Column for the base for encoding information (Required)
    --out-dir: Output folder including encoded files (Default: ./)

e.g. genofunc strain_encoder --in-fasta test/data/sequences/seqB.fasta --in-metadata test/data/metadata/metadataB.csv --encoding-column 2 --out-dir test/data/output/



### HIV Sequence Processing

Mini tools for manipulation of sequences mainly dealing with processing of raw sequences into information rich data such as annotation, feature extraction, referencing and alignment.


#### feature_extractor

Description: Extract gene(s) sequences from annotated json file based on user input gene regions. 

Options:
    --annotated-file: Input annotated json file containing all sequences (Required)
    --gene: Gene(s) regions to be extracted (Required)
    --strip-gap: Strip gaps within gene regions (Default: False)
    --filter-span: Minimum gene sequence length to be filtered (Default: 0)
    --out-prefix: Output prefix for output sequences (Default: extracted_)

e.g. genofunc feature_extractor --annotated-file test/data/sequences/annotated.json --gene pol gag --strip-gap --filter-span 0.9 --out-prefix test/data/output/extracted_ 

#### genome_annotator

Description: Annotate genomes based on closest reference sequence annotation. 

Options:
    --raw-fasta: Raw sequences with reference in name tag in fasta format (Required)
    --reference-sequence: Annotated reference sequences in json format (Required)
    --out-annotation: Output list of sequences annotated in json format (Default: referenced.json)

e.g. genofunc genome_annotator --raw-fasta test/data/sequences/referenced.fasta --reference-sequence test/data/reference/reference.json --out-annotation test/data/output/annotated.json

#### reference_matcher

Description: Map sequence to the closest reference sequence list based on mini-map2. 

Options:
    --in-fasta: Raw sequences needed to be referenced to reference list in fasta format (Required)
    --reference-sequence: Reference list in fasta format (Required)
    --out-fasta: Output list of sequences referenced (Default: referenced.fasta)

e.g. genofunc reference_matcher --in-fasta test/data/sequences/seqA.fasta --referenced-sequence test/data/reference/reference.fasta --out-fasta test/data/output/referenced.fasta

#### group_align

Description: Split the fasta file into groups of sequences set by a user threshold and align them in groups against reference. post group aligned sequences will be concatenated into a single alignment. 

Options:
    --in-dir: Input directory (Required)
    --group-size: Group size for fasta file to split by (Required)
    --reference-dir: Reference sequence directory. Reference sequence used based on matching sequence and reference file names (Required)
    --out-dir: Output folder directory (Default: ./)

e.g. genofunc group_align --in-dir test/data/sequences/group_align/ --group-size 10 --reference-dir test/data/reference/ --out-dir test/data/output/

