"""
This file is part of genofunc (https://github.com/xiaoyu518/genofunc).
Copyright 2022 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import argparse

import genofunc
import genofunc.subcommands

def main(args=None):
    parser = argparse.ArgumentParser(
        prog="genofunc",
        description="Miscellaneous genome manipulation tools",
    )

    parser.add_argument("--version", action="version", version=genofunc.__version__)
    subparsers = parser.add_subparsers(
        title="Available subcommands", help="", metavar=""
    )

    # _______________________________   common  _________________________________#
    common = argparse.ArgumentParser(add_help=False)
    
    common.add_argument(
        "-v", "--verbose", dest="verbose", action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )
    common.add_argument(
        "--log-file", dest="log_file", metavar='<filename>', required=False, default=None,
        help="Log file to use (otherwise uses stdout, or stderr if out-fasta to stdout)"
    )
    # _______________________________  concatenate_fasta  __________________________________#

    subparser_concatenate_fasta = subparsers.add_parser(
        "concatenate_fasta",
        parents=[common],
        help="Concatenate all fasta file with the same gene regions into one sequence fasta file.",
    )

    subparser_concatenate_fasta.add_argument(
        '--in-prefix', dest='in_prefix', metavar='<string>', required=True,
        help='Insert input fasta filename prefix'
    )
    subparser_concatenate_fasta.add_argument(
        '--gene', dest='gene_region', nargs='+', metavar='<list>', required=True,
        help='Insert genes to concatenate by'
    )
    subparser_concatenate_fasta.add_argument(
        '--out-dir', dest='out_dir', metavar='<foldername>', required=False, default="./",
        help='A output directory'
    )

    subparser_concatenate_fasta.set_defaults(func=genofunc.subcommands.concatenate_fasta.run)

    # _______________________________  extract_metadata  __________________________________#

    subparser_extract_metadata = subparsers.add_parser(
        "extract_metadata",
        parents=[common],
        help="Extract relevent metadata based on index fields from metadata file using fasta files. If field is empty, it will not be extracted and flagged in log file",
    )

    subparser_extract_metadata.add_argument(
        '--in-fasta', dest='in_fasta', metavar='<filename>', required=True,
        help='a FASTA file containing all sequences in metadata'
    )
    subparser_extract_metadata.add_argument(
        '--in-metadata', dest='in_metadata', metavar='<filename>', required=True,
        help='One metadata file in CSV or TSV format for extraction'
    )
    subparser_extract_metadata.add_argument(
        '--columns', dest='column_name', nargs='+', metavar='<column>', required=True,
        help='Insert column name for extraction'
    )
    subparser_extract_metadata.add_argument(
        '--rename', dest='column_rename', metavar='<list>', required=False, default=[],
        help='Insert column name for re-naming (Correspond to order of input within --c command)'
    )
    subparser_extract_metadata.add_argument(
        '--id-column', dest='id_column', metavar='<filename>', required=True,
        help='Insert column name containing sequence ID'
    )
    subparser_extract_metadata.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="filtered.fasta",
        help='Output filtered fasta file (Default: filtered.fasta)'
    )
    subparser_extract_metadata.add_argument(
        '--out-metadata', dest='out_metadata', metavar='<filename>', required=False, default="extracted_metadata.csv",
        help='Output extracted metadata file (Default: extracted_metadata.csv)'
    )

    subparser_extract_metadata.set_defaults(func=genofunc.subcommands.extract_metadata.run)

    # _______________________________  feature_extractor  __________________________________#

    subparser_feature_extractor = subparsers.add_parser(
        "feature_extractor",
        parents=[common],
        help="Extract gene sequences from referenced json file based on input gene regions.",
    )

    subparser_feature_extractor.add_argument(
        '--annotated-file', dest='in_annotation', metavar='<filename>', required=True,
        help='Annotated json file containing all sequences'
    )
    subparser_feature_extractor.add_argument(
        '--gene', dest='gene_region', nargs='+', metavar='<list>', required=True,
        help='Gene regions to be extracted'
    )
    subparser_feature_extractor.add_argument(
        '--strip-gap', dest='strip_gap', action='store_true', required=False,
        help='Strip gap bases within extracted gene regions.'
    )
    subparser_feature_extractor.add_argument(
        '--filter-span', dest='filter_span', type=float, metavar='<float>', required=False, default=0,
        help='Minimum gene sequence length to be filtered (Default: 0)'
    )
    subparser_feature_extractor.add_argument(
        '--out-prefix', dest='out_prefix', metavar='<filename>', required=False, default="extracted_",
        help='Output prefix for output sequences (Default: extracted_)'
    )

    subparser_feature_extractor.set_defaults(func=genofunc.subcommands.feature_extractor.run)

    # _______________________________  filter_fasta  __________________________________#

    subparser_filter_fasta = subparsers.add_parser(
        "filter_fasta",
        parents=[common],
        help="Filter fasta file based on minimum length.",
    )

    subparser_filter_fasta.add_argument(
        '--in-dir', dest='in_dir', metavar='<foldername>', required=True,
        help='a folder directory containing all FASTA files for filtering'
    )
    subparser_filter_fasta.add_argument(
        '--in-metadata', dest='in_metadata', metavar='<filename>', required=True,
        help='a metadata file in CSV or TSV format for filtering alongside the sequences files'
    )
    subparser_filter_fasta.add_argument(
        '--genes', dest='gene_list', nargs='+', metavar='<list>', required=True,
        help='Which genes fasta files are needed for filtering'
    )
    subparser_filter_fasta.add_argument(
        '--min-length', dest='min_length', type=float, metavar='<float>', required=True,
        help='Minimum percentage for sequences to be filtered if under this threshold'
    )
    subparser_filter_fasta.add_argument(
        '--symmetric', dest='symmetric', action='store_true', required=False,
        help='Require all gene regions to be available for the same sequence'
    )
    subparser_filter_fasta.add_argument(
        '--out-dir', dest='out_dir', metavar='<foldername>', required=False, default="./",
        help='Output directory after filtering (Default: ./)'
    )

    subparser_filter_fasta.set_defaults(func=genofunc.subcommands.filter_fasta.run)

    # _______________________________  gene_concatenator  __________________________________#

    subparser_gene_concatenator = subparsers.add_parser(
        "gene_concatenator",
        parents=[common],
        help="Concatenate fasta file based on similar sequence names of multiple genomic regions",
    )

    subparser_gene_concatenator.add_argument(
        '--in-fasta', dest='in_fasta', nargs="+", metavar='<filename>', required=True,
        help='Multiple FASTA files of sequences for concatenation'
    )
    subparser_gene_concatenator.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="concatenate.fasta",
        help='A combined single fasta file containing the combined fasta files from input'
    )

    subparser_gene_concatenator.set_defaults(func=genofunc.subcommands.gene_concatenator.run)

    # _______________________________  genome_annotator  __________________________________#

    subparser_genome_annotator = subparsers.add_parser(
        "genome_annotator",
        parents=[common],
        help="Annotate genomes based on reference sequence.",
    )

    subparser_genome_annotator.add_argument(
        '--raw-fasta', dest='raw_fasta', metavar='<filename>', required=True,
        help='One or more CSV or TSV tables of metadata'
    )
    subparser_genome_annotator.add_argument(
        '--reference-sequence', dest='reference_sequence', metavar='<filename>', required=True,
        help='Annotated reference sequences in json format'
    )
    subparser_genome_annotator.add_argument(
        '--out-annotation', dest='annotated_json', metavar='<filename>', required=False, default="referenced.json",
        help='Output list of sequences annotated in json format'
    )

    subparser_genome_annotator.set_defaults(func=genofunc.subcommands.genome_annotator.run)

    # _______________________________  group_align  __________________________________#

    subparser_group_align = subparsers.add_parser(
        "group_align",
        parents=[common],
        help="Split the fasta file into groups of sequences set by a user threshold and align them in groups.",
    )

    subparser_group_align.add_argument(
        '--in-dir', dest='in_dir', metavar='<foldername>', required=True,
        help='A directory containing FASTA files ready for group alignment'
    )
    subparser_group_align.add_argument(
        '--group-size', dest='group_size', type=int, metavar='<integer>', required=True,
        help='Group size for sequences to divide into for group alignment'
    )
    subparser_group_align.add_argument(
        '--reference-dir', dest='reference_directory', metavar='<foldername>', required=True,
        help='Group size for sequences to divide into for group alignment'
    )
    subparser_group_align.add_argument(
        '--out-dir', dest='out_dir', metavar='<foldername>', required=False, default="./", 
        help='Insert output directory (Default: ./)'
    )

    subparser_group_align.set_defaults(func=genofunc.subcommands.group_align.run)

    # _______________________________  merge  __________________________________#

    subparser_merge = subparsers.add_parser(
        "merge",
        parents=[common],
        help="Merges two or more fasta files avoiding duplicates based on matches to metadata",
    )

    subparser_merge.add_argument(
        '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=True,
        help='One or more FASTA files of sequences for merging'
    )
    subparser_merge.add_argument(
        '--in-metadata', dest='in_metadata', nargs='+', metavar='<filename>', required=True,
        help='One or more metadata in CSV or TSV format for merging'
    )
    subparser_merge.add_argument(
        '--index-field', dest='index_field', metavar='<field>', required=True,
        help='Index field containing the strain name for mergring files'
    )
    subparser_merge.add_argument(
        '--out-dir', dest='out_dir', metavar='<column>', required=False, default="./",
        help='Output directory for output fasta and metadata files'
    )

    subparser_merge.set_defaults(func=genofunc.subcommands.merge.run)

    # _______________________________  name_splitter  __________________________________#

    subparser_name_splitter = subparsers.add_parser(
        "name_splitter",
        parents=[common],
        help="Split the sequence name into metadata based on piping character."
    )

    subparser_name_splitter.add_argument(
        '--in-fasta', dest='in_fasta', metavar='<filename>', required=True,
        help='Input fasta file'
    )
    subparser_name_splitter.add_argument(
        '--pipe', dest='piping_character', metavar='<character>', required=True,
        help='Input character the fasta sequence name is split based on'
    )
    subparser_name_splitter.add_argument(
        '--header', dest='header', nargs='+', metavar='<list>', required=False, default="",
        help='Header for the output metadata table.'
    )
    subparser_name_splitter.add_argument(
        '--out-metadata', dest='out_metadata', metavar='<filename>', required=False, default="metadata.csv",
        help='Output metadata file (Default: metadata.csv)'
    )

    subparser_name_splitter.set_defaults(func=genofunc.subcommands.name_splitter.run)


    # _______________________________  reference_matcher  __________________________________#

    subparser_reference_matcher = subparsers.add_parser(
        "reference_matcher",
        parents=[common],
        help="Map sequence to reference sequence list."
    )

    subparser_reference_matcher.add_argument(
        '--in-fasta', dest='in_fasta', metavar='<filename>', required=True,
        help='Raw sequences needed to be referenced to reference list in fasta format'
    )
    subparser_reference_matcher.add_argument(
        '--reference-sequence', dest='reference_sequence', metavar='<filename>', required=True,
        help='Reference list in fasta format'
    )
    subparser_reference_matcher.add_argument(
        '--out_fasta', dest='out_fasta', default="referenced.fasta", metavar='<filename>', required=False,
        help='Output list of referenced sequences (Default:referenced.fasta)'
    )

    subparser_reference_matcher.set_defaults(func=genofunc.subcommands.reference_matcher.run)

    # _______________________________  rename_fasta  __________________________________#

    subparser_rename_fasta = subparsers.add_parser(
        "rename_fasta",
        parents=[common],
        help="Renaming fasta sequence names based on character splits"
    )

    subparser_rename_fasta.add_argument(
        '--in-fasta', dest='in_fasta', metavar='<filename>', required=True,
        help='Input sequence fasta file'
    )

    subparser_rename_fasta.add_argument(
        '--pipe', dest='piping_character', metavar='<character>', required=True,
        help='Input character the fasta sequence name is split based on'
    )

    subparser_rename_fasta.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="cleaned_sequences.fasta",
        help='Output fasta file (Default: cleaned_sequences.fasta)'
    )

    subparser_rename_fasta.set_defaults(func=genofunc.subcommands.rename_fasta.run)

    # _______________________________  strain_encoder  __________________________________#

    subparser_strain_encoder = subparsers.add_parser(
        "strain_encoder",
        parents=[common],
        help="Encoded strain id into non-defining ids"
    )

    subparser_strain_encoder.add_argument(
        '--in-fasta', dest='in_fasta', metavar='<filename>', required=True,
        help='Input folder containing all fasta files for encoding'
    )

    subparser_strain_encoder.add_argument(
        '--in-metadata', dest='in_metadata', metavar='<filename>', required=True,
        help='Input metadata in CSV or TSV format'
    )
    subparser_strain_encoder.add_argument(
        '--encoding-column', dest='encoding_column', type=int, metavar='<integer>', required=True,
        help='Column number for the base of encoding information'
    )

    subparser_strain_encoder.add_argument(
        '--out-dir', dest='out_dir', metavar='<foldername>', required=False, default="./",
        help='Output folder including encoded files'
    )

    subparser_strain_encoder.set_defaults(func=genofunc.subcommands.strain_encoder.run)

    # ___________________________________________________________________________________#

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
