FASTQGEN
========

A simple tool to generate random paired-end FASTQ files for testing and
development purposes.


DESCRIPTION
-----------

fastqgen creates synthetic paired-end sequencing reads in FASTQ format.
The tool generates:

- Random DNA sequences (A, T, C, G)
- Reverse complement mate pairs
- Phred quality scores (Q0-Q40, ASCII 33-73)
- Properly formatted FASTQ output files


INSTALLATION
------------

**Prerequisites:**

You need Rust and Cargo installed. If you use conda/mamba:

    mamba install -c conda-forge rust
    # or
    conda install -c conda-forge rust

Alternatively, install Rust from https://rustup.rs

**From crates.io**

    cargo install fastqgen

**From Github**

    cargo install --path .

This will compile the binary and install it to ~/.cargo/bin/
Make sure ~/.cargo/bin is in your PATH.


USAGE
-----

Basic usage:

    fastqgen -n 1000 -l 150 -o my_reads

Options:

    -o, --outfile <NAME>    Output file prefix [default: synthetic_reads]
    -n <NUMBER>             Number of read pairs to generate [default: 1000]
    -l <LENGTH>             Read length in base pairs [default: 150]
    -h, --help              Print help
    -V, --version           Print version


OUTPUT
------

The tool generates two files:

    <outfile>_R1.fastq      Forward reads
    <outfile>_R2.fastq      Reverse reads (reverse complement of R1)

Each FASTQ record contains:
- Header line with read ID and pair indicator (/1 or /2)
- Sequence line
- Plus line separator
- Quality score line


EXAMPLES
--------

Generate 5000 read pairs of 100bp length:

    fastqgen -n 5000 -l 100 -o test_data

Generate default dataset (1000 reads, 150bp):

    fastqgen


LICENSE
-------

MIT


AUTHOR
------

Daniel A. Sprague