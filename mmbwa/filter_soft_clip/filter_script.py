import sys
import pysam
import typer

app = typer.Typer()

def total_soft_clip_length(read):
    cig = read.cigartuples
    if cig is None:
        return 0
    return sum(length for op, length in cig if op == 4)

def write_fastq(read):
    name = read.query_name
    seq = read.query_sequence
    qual = read.query_qualities
    if seq is None or qual is None:
        typer.echo(f"Warning: Skipped read {name} due to missing sequence or qualities", err=True)
        return
    qual_str = "".join(chr(q + 33) for q in qual)
    sys.stdout.write(f"@{name}\n{seq}\n+\n{qual_str}\n")
    sys.stdout.flush()

def filter_alns(
    bam_path=None,
    output_bam_path=None,
    output_temp=None,
    threshold_fraction=0.10,
    unmapped=False,
    keep_temp=False
):
    # Read BAM from file or stdin
    if bam_path is None:
        bam = pysam.AlignmentFile("-", "r")  # SAM from stdin
    else:
        if bam_path.endswith(".bam"):
            bam = pysam.AlignmentFile(bam_path, "rb") # read BAM from pre-aligned file
        else:
            bam = pysam.AlignmentFile(bam_path, "r")  #read SAM from pre-aligned file

    # Prepare output BAM for filtered reads (below threshold)
    output_bam = pysam.AlignmentFile(output_bam_path, "wb", header=bam.header)
    # If temporary files are kept write it to the same output path
    temp_bam = None
    if keep_temp:
        temp_bam = pysam.AlignmentFile(output_temp, "wb", header=bam.header)

    for read in bam:
        if read.is_secondary or read.is_supplementary: # process only primary alignments
            continue
        if read.query_length is None: # edge case
            continue
        if not unmapped and read.is_unmapped: # filter unmapped reads if --unmapped
            continue

        total_soft_clip = total_soft_clip_length(read)

        if total_soft_clip > threshold_fraction * read.query_length:
            # write heavily soft clipped reads to stdout or to temp bam file
            if keep_temp:
                temp_bam.write(read)
                write_fastq(read)
            else:
                write_fastq(read)
        else:
            output_bam.write(read)

    if temp_bam:
        temp_bam.close()

    bam.close()
    output_bam.close()

@app.command()
def main(
    input_bam: str = typer.Argument(
        None,
        help="Input BAM file path. If omitted, reads BAM from stdin.",
    ),
    output_bam: str = typer.Option(
        ...,
        "-o",
        "--output",
        help="Output BAM file path for reads below soft clipping threshold",
    ),
    output_temp: str = typer.Option(None, "--output-temp", help="Temporary BAM file for realignment"),
    threshold: float = typer.Option(
        0.10,
        help="Soft clipping fraction threshold to filter_soft_clip reads",
    ),
    unmapped: bool = typer.Option(
        False,
        help="Align unmapped reads with bwa mem", ),
    keep_temp: bool = typer.Option(
            False,
            help="Keep temporary files", ),
):
    filter_alns(
        bam_path=input_bam,
        output_bam_path=output_bam,
        output_temp=output_temp,
        threshold_fraction=threshold,
        unmapped = unmapped,
        keep_temp = keep_temp,
    )

if __name__ == "__main__":
    app()
