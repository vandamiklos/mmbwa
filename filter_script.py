import sys
import pysam
import typer
from collections import defaultdict

app = typer.Typer()

def soft_clip_length(read):
    if read.cigartuples is None:
        return 0
    cigar = read.cigartuples
    soft_clip_start = cigar[0][1] if cigar[0][0] == 4 else 0
    soft_clip_end = cigar[-1][1] if cigar[-1][0] == 4 else 0
    return soft_clip_start + soft_clip_end

def write_fastq(read):
    name = read.query_name
    seq = read.query_sequence
    qual = read.query_qualities
    if seq is None or qual is None:
        return
    qual_str = "".join(chr(q + 33) for q in qual)
    sys.stdout.write(f"@{name}\n{seq}\n+\n{qual_str}\n")

def filter_reads_by_total_soft_clipping(
    bam_path=None,
    output_bam_path=None,
    threshold_fraction=0.1,
    write_fastq_flag=False
):
    # Read BAM from file or stdin
    if bam_path is None:
        bam = pysam.AlignmentFile("-", "rb")  # BAM from stdin
    else:
        bam = pysam.AlignmentFile(bam_path, "rb")

    # Prepare output BAM for filtered reads (below threshold)
    output_bam = pysam.AlignmentFile(output_bam_path, "wb", header=bam.header)

    read_segments = defaultdict(list)
    for read in bam:
        if read.is_unmapped:
            continue
        read_segments[read.query_name].append(read)
    bam.close()

    count_written = 0
    for qname, segments in read_segments.items():
        read_length = sum(r.query_length for r in segments if r.query_length)
        if read_length is None:
            continue

        total_soft_clip = sum(soft_clip_length(r) for r in segments)

        if total_soft_clip > threshold_fraction * read_length:
            # Write heavily soft clipped reads to FASTQ if requested
            if write_fastq_flag:
                for r in segments:
                    write_fastq(r)
        else:
            # Write below threshold reads to output BAM
            for r in segments:
                output_bam.write(r)
            count_written += 1

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
    threshold: float = typer.Option(
        0.10,
        help="Soft clipping fraction threshold to filter reads",
    ),
    write_fastq: bool = typer.Option(
        False,
        help="Write reads above threshold to FASTQ on stdout",
    ),
):
    filter_reads_by_total_soft_clipping(
        bam_path=input_bam,
        output_bam_path=output_bam,
        threshold_fraction=threshold,
        write_fastq_flag=write_fastq,
    )

if __name__ == "__main__":
    app()
