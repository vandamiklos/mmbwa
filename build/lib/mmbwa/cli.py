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
    input_fastq: Path = typer.Option(None, "--input-fq", help="Input FASTQ file"),
    input_aln: Path = typer.Option(None, "--input-aln", help="Input SAM/BAM file if already aligned by minimap2"),
    threads: int = typer.Option(1, "--threads", "-t"),
    mm_args: str = typer.Option("-ax map-ont", "--mm-args"),
    bwa_args: str = typer.Option("", "--bwa-args"),
    output: Path = typer.Option(..., "--output", "-o"),
    filter_script: str = typer.Option("filter_script.py"),
    keep_temp: bool = typer.Option(False, "--keep-temp"),
    threshold: float = typer.Option(0.10, "--threshold"),
    sort: bool = typer.Option(False, "--sort"),
    index: bool = typer.Option(False, "--index"),
    unmapped: bool = typer.Option(False, "--unmapped")
):
    output.mkdir(parents=True, exist_ok=True)
    realigned_bam = output / "realigned.bam"
    sorted_bam = output / "final.sorted.bam"
    filtered_bam = output / "filtered.bam"

    if input_aln:
        # Use prealigned BAM/SAM as input to filter_script and bwa mem
        # Filter reads by soft clipping and output filtered BAM and soft-clipped FASTQ for bwa mem
        filter_cmd = f"python {filter_script} {input_aln} --output {filtered_bam} --output-temp /dev/null --threshold {threshold}"
        if unmapped:
            filter_cmd += " --unmapped"
        if keep_temp:
            # If you want to keep temp files, define output_temp path
            softclipped_bam = output / "softclipped.bam"
            filter_cmd += f" --keep-temp --output-temp {softclipped_bam}"
            bwa_input = softclipped_bam
            bwa_cmd = f"bwa mem {bwa_args} -t {threads} {ref} {bwa_input} | samtools view -b -o {realigned_bam} -"
            # Run filter first, then bwa mem on softclipped.bam
            run_cmd(filter_cmd)
            run_cmd(bwa_cmd)
        else:
            # Pipe filter stdout (FASTQ) to bwa mem
            bwa_cmd = f"bwa mem {bwa_args} -t {threads} {ref} - | samtools view -b -o {realigned_bam} -"
            full_cmd = f"{filter_cmd} | {bwa_cmd}"
            run_cmd(full_cmd)

    else:
        # minimap2 | filter_script | bwa mem | samtools
        if input_fastq is None:
            typer.echo("Error: Either --input-fq or --input-aln must be provided.", err=True)
            raise typer.Exit(code=1)

        if keep_temp:
            mm_sam = output / "minimap2.sam"
            softclipped_bam = output / "softclipped.bam"
            mm_cmd = f"minimap2 {mm_args} -t {threads} {ref} {input_fastq} > {mm_sam}"
            run_cmd(mm_cmd)

            filter_cmd = f"python {filter_script} {mm_sam} --output {filtered_bam} --output-temp {softclipped_bam} --threshold {threshold}"
            if unmapped:
                filter_cmd += " --unmapped"
            if keep_temp:
                filter_cmd += " --keep-temp"

            bwa_cmd = f"bwa mem {bwa_args} -t {threads} {ref} {softclipped_bam} | samtools view -b -o {realigned_bam}"

            full_cmd = f"{filter_cmd} | {bwa_cmd}"
            run_cmd(full_cmd)
        else:
            mm_cmd = f"minimap2 {mm_args} -t {threads} {ref} {input_fastq}"
            filter_cmd = f"python {filter_script} - --output {filtered_bam} --output-temp /dev/null --threshold {threshold}"
            if unmapped:
                filter_cmd += " --unmapped"

            bwa_cmd = f"bwa mem {bwa_args} -t {threads} {ref} - | samtools view -b -o {realigned_bam}"

            full_cmd = f"{mm_cmd} | {filter_cmd} | {bwa_cmd}"
            run_cmd(full_cmd)

        # Optional sort/index

    final_bam = output / "final.bam"
    run_cmd(f"samtools merge -@ {threads} {final_bam} {filtered_bam} {realigned_bam}")

    # Optional sort/index
    if sort:
        run_cmd(f"samtools sort -@ {threads} -o {sorted_bam} {final_bam}")
        final_bam = sorted_bam
    if index:
        run_cmd(f"samtools index {final_bam}")

    typer.echo(f"Output BAM: {final_bam}")

    if not keep_temp:
        for f in [realigned_bam, filtered_bam, softclipped_bam, mm_sam]:
            if f.exists():
                f.unlink()
                typer.echo(f"Deleted temporary file: {f}")

def main():
    app()

if __name__ == "__main__":
    main()
