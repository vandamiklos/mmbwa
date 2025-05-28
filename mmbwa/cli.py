# mmbwa/cli.py

import typer
import subprocess
from pathlib import Path

app = typer.Typer(help="Realign soft-clipped reads with bwa mem after minimap2.")

def run_cmd(command: str):
    typer.echo(f"Running: {command}")
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        typer.echo(f"Command failed: {command}", err=True)
        raise typer.Exit(code=result.returncode)

@app.command()
def realign(
    ref: Path = typer.Argument(..., help="Reference genome FASTA file"),
    input_fastq: Path = typer.Argument(..., help="Input FASTQ file"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
    mm_args: str = typer.Option("-ax map-ont", "--mm-args", help="Extra minimap2 arguments"),
    bwa_args: str = typer.Option("", "--bwa-args", help="Extra bwa mem arguments"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    filter_script: str = typer.Option("filter_script.py", help="Script to filter soft-clipped reads"),
    keep_intermediate: bool = typer.Option(False, "--keep-intermediate", help="Keep intermediate BAM files"),
    threshold: float = typer.Option(0.10, "--threshold", help="Soft clipping fraction threshold for filtering reads"),
):
    """
    Re-align soft-clipped reads detected by minimap2 using bwa mem.
    Writes only final.sorted.bam and its index unless --keep-intermediate is specified.
    """
    output.mkdir(parents=True, exist_ok=True)
    filtered_bam = output / "filtered.bam"
    realigned_bam = output / "realigned.bam"
    final_bam = output / "final.bam"
    sorted_bam = output / "final.sorted.bam"

    # Step 1: Run minimap2, convert to BAM, pipe into filter_script.py
    mm_cmd = f"minimap2 {mm_args} -t {threads} {ref} {input_fastq} | samtools view -b -"
    filter_cmd = f"python {filter_script} -o {filtered_bam} --threshold {threshold} --write-fastq"

    # Step 2: Pipe FASTQ output of filter_script.py into bwa mem, convert to BAM
    bwa_mem_cmd = f"bwa mem {bwa_args} -t {threads} {ref} - | samtools view -b - > {realigned_bam}"

    # Combine all commands with pipes
    full_cmd = f"{mm_cmd} | {filter_cmd} | {bwa_mem_cmd}"
    run_cmd(full_cmd)

    # Step 3: Merge filtered BAM and realigned BAM
    merge_cmd = f"samtools merge -@ {threads} -o {final_bam} {filtered_bam} {realigned_bam}"
    run_cmd(merge_cmd)

    # Step 4: Sort the merged BAM
    sort_cmd = f"samtools sort -@ {threads} -o {sorted_bam} {final_bam}"
    run_cmd(sort_cmd)

    # Step 5: Index the sorted BAM
    index_cmd = f"samtools index {sorted_bam}"
    run_cmd(index_cmd)

    typer.echo(f"Sorted and indexed BAM: {sorted_bam}")

    # Cleanup intermediate files if not keeping them
    if not keep_intermediate:
        for f in [filtered_bam, realigned_bam, final_bam]:
            if f.exists():
                f.unlink()
                typer.echo(f"Deleted intermediate file: {f}")

def main():
    app()
