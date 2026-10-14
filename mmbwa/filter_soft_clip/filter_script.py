import sys
import pysam
import typer
import numpy as np
from sortedintersect import IntervalSet
from collections import defaultdict, namedtuple
import pandas as pd

app = typer.Typer()


def soft_clip_length(read):
    if read.cigartuples is None or len(read.cigartuples) == 0:
        return 0
    softclip = 0
    if read.cigartuples[0][0] == pysam.CSOFT_CLIP:
        softclip += read.cigartuples[0][1]
    if read.cigartuples[-1][0] == pysam.CSOFT_CLIP:
        softclip += read.cigartuples[-1][1]
    return softclip


def write_fastq(read):
    name = read.query_name
    seq = read.query_sequence
    qual = read.query_qualities
    if seq is None or qual is None:
        typer.echo(f"Warning: Skipped read {name} due to missing sequence or qualities", err=True)
        return
    qual_str = pysam.array_to_qualitystring(qual)
    sys.stdout.write(f"@{name}\n{seq}\n+\n{qual_str}\n")


IntervalItem = namedtuple('interval_item',
                          ['chrom', 'start', 'end'])


def prepare_data_from_bed(bed_df, col):
    # need to make sure rend > rstart for sortedintersect, and intervals are sorted
    bed_df['start'] = np.minimum(bed_df.iloc[:, 1], bed_df.iloc[:, 2])
    bed_df['end'] = np.maximum(bed_df.iloc[:, 1], bed_df.iloc[:, 2])
    bed_df['chrom'] = bed_df.iloc[:, 0]
    bed_df = bed_df.sort_values('start')
    data = []
    for i in zip(*(bed_df[c] for c in col)):
        data.append(IntervalItem(*i))
    return data


def build_interval_trees(data):
    interval_tree = defaultdict(lambda: IntervalSet(with_data=True))
    for itv in data:
        interval_tree[itv.chrom].add(itv.start, itv.end, itv)
    return interval_tree


def filter_long_soft_clip(
    bam_path=None,
    output_bam_path=None,
    output_temp=None,
    threshold=100,
    unmapped=False,
    keep_temp=False,
    regions_bed=None,
    len_filter=None,
    mapq_filter=None
):
    # Read BAM from file or stdin
    if bam_path is None:
        bam = pysam.AlignmentFile("-", "r")  # SAM from stdin
    else:
        if bam_path.endswith(".bam"):
            bam = pysam.AlignmentFile(bam_path, "rb")  # read BAM from pre-aligned file
        else:
            bam = pysam.AlignmentFile(bam_path, "r")  # read SAM from pre-aligned file

    # Prepare output BAM for filtered reads (below threshold)
    output_bam = pysam.AlignmentFile(output_bam_path, "wb", header=bam.header)
    # If temporary files are kept write it to the same output path
    temp_bam = None
    if keep_temp:
        temp_bam = pysam.AlignmentFile(output_temp, "wb", header=bam.header)

    interval_tree = {}

    if regions_bed is not None:
        bed = pd.read_csv(regions_bed,sep="\t",header=None,comment="#",usecols=[0, 1, 2],)
        bed.columns = ["chrom", "start", "end"]
        data = prepare_data_from_bed(
            bed,["chrom", "start", "end"],)
    interval_tree = build_interval_trees(data)

    for read in bam:
        if read.is_secondary or read.is_supplementary:  # process only primary alignments
            continue
        if read.query_length is None:  # edge case
            continue
        if read.is_unmapped:
            if unmapped:
                if keep_temp:
                    temp_bam.write(read)
                write_fastq(read)
            continue
        if read.query_sequence is None:
            output_bam.write(read)
            continue
        if len(read.query_sequence) > len_filter:
            output_bam.write(read)
            continue
        if read.mapq > mapq_filter:
            output_bam.write(read)
            continue


        overlap_interval = interval_tree[read.reference_name].search_interval(read.reference_start, read.reference_end)
        if overlap_interval:
            output_bam.write(read)
            continue

        if soft_clip_length(read) > threshold:
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
    input_bam: str = typer.Argument(None, help="Input BAM file path. If omitted, reads BAM from stdin.",),
    output_bam: str = typer.Option(..., "-o", "--output", help="Output BAM file path for reads below soft "
                                                               "clipping threshold",),
    output_temp: str = typer.Option(None, "--output-temp", help="Temporary BAM file for realignment"),
    threshold: int = typer.Option(100, help="Soft clipping length threshold to filter_soft_clip reads", ),
    unmapped: bool = typer.Option(False, help="Align unmapped reads with bwa mem", ),
    keep_temp: bool = typer.Option(False, help="Keep temporary files", ),
    regions_bed: str = typer.Option(None, "--regions_bed", "-rb", help="Genomic regions in BED format"),
    len_filter: int = typer.Option(100_000, help="Don't realign reads that are longer than this threshold "),
    mapq_filter: int = typer.Option(30, help="Don't realign if above this mapq"),
):
    filter_long_soft_clip(
        bam_path=input_bam,
        output_bam_path=output_bam,
        output_temp=output_temp,
        threshold=threshold,
        unmapped=unmapped,
        keep_temp=keep_temp,
        regions_bed=regions_bed,
        len_filter=len_filter,
        mapq_filter=mapq_filter
    )


if __name__ == "__main__":
    app()
