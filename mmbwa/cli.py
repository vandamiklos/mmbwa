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
    input_fq: Path = typer.Option(None, "--input-fq", help="Input FASTQ file"),
    input_aln: Path = typer.Option(None, "--input-aln", help="Input SAM/BAM file if already aligned by minimap2"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
    mm_args: str = typer.Option("-ax map-ont", "--mm-args", help="Arguments passed to minimap2"),
    bwa_args: str = typer.Option("", "--bwa-args", help="Arguments passed to bwa mem"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    keep_temp: bool = typer.Option(False, "--keep-temp", help="Keep temporary files"),
    threshold: float = typer.Option(100, "--threshold", help='Length of soft-clips'),
    sort: bool = typer.Option(False, "--sort", help='Sort final bam output'),
    index: bool = typer.Option(False, "--index", help='Generate index file from the sorted final bam output'),
    unmapped: bool = typer.Option(False, "--unmapped", help='Re-align unmapped reads with bwa mem'),
    log: bool = typer.Option(False, "--log", help='Output log file'),
    no_merge: bool = typer.Option(False, "--no-merge", help='Leave minimap2 and bwa-mem alignment files separately')
):
    output.mkdir(parents=True, exist_ok=True)
    realigned_bam = output / "realigned.bam"
    sorted_bam = output / "final.sorted.bam"
    filtered_bam = output / "filtered.bam"
    softclipped_bam = output / "softclipped.bam" if keep_temp else None
    mm_sam = output / "minimap2.sam" if keep_temp else None
    log_file = output / "log.txt" if log else None

    if log:
        log = open(log_file, "a")

# if pre-aligned minimap2 SAM/BAM is available: BAM/SAM -> filter_soft_clip | bwa-mem
    if input_aln:

        filter_args = [
            "python", "-m", "mmbwa.filter_soft_clip.filter_script",
            str(input_aln),
            "--output", str(filtered_bam),
            "--threshold", str(threshold)]

        if keep_temp:
            filter_args += ["--output-temp", str(softclipped_bam), "--keep-temp"]
        else:
            filter_args += ["--output-temp", "/dev/null"]

        if unmapped:
            filter_args.append("--unmapped")

        bwa_args_list = ["bwa", "mem"] + bwa_args.strip().split() + ["-t", str(threads), str(ref), "-"]
        samtools_args = ["samtools", "view", "-b", "-o", str(realigned_bam)]

        filter_proc = subprocess.Popen(filter_args, stdout=subprocess.PIPE, stderr=log, text=True)
        bwa_proc = subprocess.Popen(bwa_args_list, stdin=filter_proc.stdout, stdout=subprocess.PIPE, stderr=log,
                                    text=True)

        # close filter_proc.stdout, so it gets a SIGPIPE if bwa_proc exits early
        filter_proc.stdout.close()

        with open(realigned_bam, "wb") as bam_out:
            subprocess.run(samtools_args, stdin=bwa_proc.stdout, stdout=bam_out, stderr=log, check=True)

    # FASTQ -> minimap2 | filter_soft_clip | bwa-mem
    else:
        if input_fq is None:
            typer.echo("Error: Either --input-fq or --input-aln must be provided.", err=True)
            raise typer.Exit(code=1)

        if keep_temp:
            # run minimap2 and tee its output
            mm_proc = subprocess.Popen(
                ["minimap2"] + mm_args.strip().split() + ["-t", str(threads), str(ref), str(input_fq)],
                stdout=subprocess.PIPE, text=True, stderr=log)

            # tee minimap2 output to a file and pass to filter
            mm_tee_proc = subprocess.Popen(
                ["tee", str(mm_sam)],
                stdin=mm_proc.stdout,
                stdout=subprocess.PIPE, text=True, stderr=log)

            mm_proc.stdout.close()

            # run filter script and tee its output
            filter_args = [
                "python", "-m", "mmbwa.filter_soft_clip.filter_script", "-",
                "--output", str(filtered_bam),
                "--output-temp", str(softclipped_bam),
                "--threshold", str(threshold),
                "--keep-temp"]

            if unmapped:
                filter_args.append("--unmapped")

            filter_proc = subprocess.Popen(
                filter_args,
                stdin=mm_tee_proc.stdout,
                stdout=subprocess.PIPE, text=True, stderr=log)

            mm_tee_proc.stdout.close()

            # BWA MEM on filtered output
            bwa_proc = subprocess.Popen(
                ["bwa", "mem"] + bwa_args.strip().split() + ["-t", str(threads), str(ref), "-"],
                stdin=filter_proc.stdout,
                stdout=subprocess.PIPE, text=True, stderr=log)

            filter_proc.stdout.close()

            with open(realigned_bam, "wb") as bam_out:
                subprocess.run(["samtools", "view", "-b"], stdin=bwa_proc.stdout, stdout=bam_out, text=True, stderr=log)

            bwa_proc.stdout.close()

        else:
            # minimap2 subprocess
            mm_proc = subprocess.Popen(
                ["minimap2"] + mm_args.strip().split() + ["-t", str(threads), str(ref), str(input_fq)],
                stdout=subprocess.PIPE, text=True, stderr=log)

            # filter_soft_clip subprocess
            filter_args = [
                "python", "-m", "mmbwa.filter_soft_clip.filter_script", "-",
                "--output", str(filtered_bam),
                "--output-temp", "/dev/null",
                "--threshold", str(threshold)]

            if unmapped:
                filter_args.append("--unmapped")

            filter_proc = subprocess.Popen(filter_args, stdin=mm_proc.stdout, stdout=subprocess.PIPE, text=True,
                                           stderr=log)

            mm_proc.stdout.close()

            # bwa mem subprocess
            bwa_proc = subprocess.Popen(
                ["bwa", "mem"] + bwa_args.strip().split() + ["-t", str(threads), str(ref), "-"],
                stdin=filter_proc.stdout,
                stdout=subprocess.PIPE
            )

            with open(realigned_bam, "wb") as bam_out:
                subprocess.run(["samtools", "view", "-b"], stdin=bwa_proc.stdout, stdout=bam_out, text=True, stderr=log)

            filter_proc.stdout.close()

    if not no_merge:
        final_bam = output / "final.bam"
        typer.echo(f"Merging {filtered_bam} and {realigned_bam}")
        run_cmd(f"samtools merge -@ {threads} {final_bam} {filtered_bam} {realigned_bam}")

        # Optional sort/index
        if sort:
            typer.echo(f"Sorting {final_bam}")
            run_cmd(f"samtools sort -@ {threads} -o {sorted_bam} {final_bam}")
            final_bam = sorted_bam
        if index:
            typer.echo(f"Indexing {final_bam}")
            run_cmd(f"samtools index {final_bam}")

        typer.echo(f"Output BAM: {final_bam}")

    if not keep_temp:
        for f in [realigned_bam, filtered_bam, softclipped_bam, mm_sam]:
            if f and f.exists():
                f.unlink()
                typer.echo(f"Deleted temporary file: {f}")

    if log:
        log.close()


def main():
    app()


if __name__ == "__main__":
    main()
