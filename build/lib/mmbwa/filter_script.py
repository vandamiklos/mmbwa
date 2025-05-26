import sys
import pysam
import argparse

def is_soft_clipped(cigar):
    return any(op in [4, 5] and length > 20 for (op, length) in cigar)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", required=True, help="Output BAM file for non-soft-clipped reads")
    args = parser.parse_args()

    # input from stdin (SAM), output BAM
    in_sam = pysam.AlignmentFile("-", "r")  # SAM from minimap2
    out_bam = pysam.AlignmentFile(args.output, "wb", header=in_sam.header)

    for read in in_sam.fetch(until_eof=True):
        if read.cigartuples and is_soft_clipped(read.cigartuples):
            sys.stdout.write(read.to_string() + "\n")  # Write SAM line to stdout
        else:
            out_bam.write(read)

    out_bam.close()

if __name__ == "__main__":
    main()
