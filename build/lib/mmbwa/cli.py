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
    input_fq: Path = typer.Option(None, "--input-fq", help="Input FASTQ file - no prior alignment available"),
    input_aln: Path = typer.Option(None, "--input-aln", help="Input SAM/BAM file if already aligned by minimap2, "
                                                             "must be sorted and indexed"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
    mm_args: str = typer.Option("-ax map-ont", "--mm-args", help="Arguments passed to minimap2"),
    bwa_args: str = typer.Option("", "--bwa-args", help="Arguments passed to bwa mem"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    regions_bed: Path = typer.Option(None, "--regions_bed", "-rb", help="Genomic regions to exclude from re-alignment "
                                                                        "in BED format"),
    regions_overlap: float = typer.Option(0.7, "--regions_overlap", "-ro", help="Min overlap with regions"),
    keep_temp: bool = typer.Option(False, "--keep-temp", help="Keep temporary files"),
    threshold: int = typer.Option(100, "--threshold", help='Threshold for soft-clip length, above this the read will '
                                                             'be filtered'),
    sort: bool = typer.Option(False, "--sort", help='Sort final bam output using samtools sort'),
    index: bool = typer.Option(False, "--index", help='Generate index file from the sorted final bam output, using '
                                                      'samtools index'),
    unmapped: bool = typer.Option(False, "--unmapped", help='Re-align unmapped reads with bwa mem'),
    log: bool = typer.Option(False, "--log", help='Output log file'),
    no_merge: bool = typer.Option(False, "--no-merge", help='Leave minimap2 and bwa mem alignment files separately')
):
    output.mkdir(parents=True, exist_ok=True)
    realigned_bam = output / "realigned.bam"
    sorted_bam = output / "final.sorted.bam"
    filtered_bam = output / "filtered.bam"
    softclipped_bam = output / "softclipped.bam" if keep_temp else None
    mm_sam = output / "minimap2.sam" if keep_temp else None
    log_file = output / "log.txt" if log else None

    log_fh = open(log_file, "a", buffering=1) if log else None

    bwa_args_list = ["bwa", "mem"] + bwa_args.strip().split() + ["-t", str(threads), str(ref), "-"]

    # if pre-aligned minimap2 SAM/BAM is available: BAM/SAM -> filter_soft_clip | bwa-mem
    if input_aln:
        filter_args = [
            "python", "-m", "mmbwa.filter_soft_clip.filter_script",
            str(input_aln),
            "--output", str(filtered_bam),
            "--threshold", str(threshold),
            "--regions_bed", str(regions_bed)]

#        if regions:
#            not_intersect_bam = output / "not_in_regions.bam"
#            intersect_bam = output / "within_regions.bam"
#            tmp_primary_bam = output / "tmp_primary.bam"
#            # save primary alignments that are not in specific regions
#            cmd = (
#                f"samtools view -h -b -F 0x900 {input_aln} > {tmp_primary_bam} ;"
#                f"bedtools intersect -abam {tmp_primary_bam} -b {regions_bed} -v -f {regions_overlap}"
#                f"| "
#                f"samtools sort -o {not_intersect_bam} -"
#            )
#
#            subprocess.run(
#                cmd,
#                shell=True,
#                stderr=log_fh,
#                check=True
#            )
#
#            cmd2 = (
#                f"bedtools intersect -abam {tmp_primary_bam} -b {regions_bed} -wa -f {regions_overlap} "
#                f"| "
#                f"samtools sort -o {intersect_bam} -"
#            )
#
#            subprocess.run(
#                cmd2,
#                shell=True,
#                stderr=log_fh,
#                check=True
#            )
#
#            # change input file
#            filter_args = [
#                "python", "-m", "mmbwa.filter_soft_clip.filter_script",
#                str(not_intersect_bam),
#                "--output", str(filtered_bam),
#                "--threshold", str(threshold)]

        if keep_temp:
            filter_args += ["--output-temp", str(softclipped_bam), "--keep-temp"]
        if not keep_temp:
            filter_args += ["--output-temp", "/dev/null"]
        if unmapped:
            filter_args.append("--unmapped")


        filter_proc = subprocess.Popen(filter_args, stdout=subprocess.PIPE, stderr=log_fh)
        bwa_proc = subprocess.Popen(bwa_args_list, stdin=filter_proc.stdout, stdout=subprocess.PIPE,
                                    stderr=log_fh)
        # close filter_proc.stdout, so it gets a SIGPIPE if bwa_proc exits early
        filter_proc.stdout.close()

        subprocess.run(["samtools", "view", "-b", "-o", str(realigned_bam)],
                           stdin=bwa_proc.stdout, stderr=log_fh, check=True)

        # ensure upstream processes are reaped and check exit codes
        bwa_rc = bwa_proc.wait()
        filter_rc = filter_proc.wait()
        if bwa_rc != 0:
            raise typer.Exit(code=bwa_rc)
        if filter_rc != 0:
            raise typer.Exit(code=filter_rc)

    # FASTQ -> minimap2 | filter_soft_clip | bwa-mem
    else:
        if input_fq is None:
            typer.echo("Error: Either --input-fq or --input-aln must be provided.", err=True)
            raise typer.Exit(code=1)

        minimap2_args_list = ["minimap2"] + mm_args.strip().split() + ["-t", str(threads), str(ref), str(input_fq)]

        if keep_temp:
            # run minimap2 and tee its output
            mm_proc = subprocess.Popen(minimap2_args_list, stdout=subprocess.PIPE, stderr=log_fh)
            # tee minimap2 output to a file and pass to filter
            mm_tee_proc = subprocess.Popen(["tee", str(mm_sam)], stdin=mm_proc.stdout, stdout=subprocess.PIPE,
                                           stderr=log_fh)
            mm_proc.stdout.close()

            # run filter script and tee its output
            filter_args = [
                "python", "-m", "mmbwa.filter_soft_clip.filter_script", "-",
                "--output", str(filtered_bam),
                "--output-temp", str(softclipped_bam),
                "--threshold", str(threshold),
                "--keep-temp",
                "--regions_bed", str(regions_bed)]

            if unmapped:
                filter_args.append("--unmapped")

            filter_proc = subprocess.Popen(filter_args, stdin=mm_tee_proc.stdout, stdout=subprocess.PIPE, stderr=log_fh)
            mm_tee_proc.stdout.close()

            # BWA MEM on filtered output
            bwa_proc = subprocess.Popen(bwa_args_list, stdin=filter_proc.stdout, stdout=subprocess.PIPE, stderr=log_fh)
            filter_proc.stdout.close()

            with open(realigned_bam, "wb") as bam_out:
                subprocess.run(["samtools", "view", "-b"], stdin=bwa_proc.stdout, stdout=bam_out,
                               stderr=log_fh)

            # ensure upstream processes are reaped and check exit codes
            bwa_rc = bwa_proc.wait()
            filter_rc = filter_proc.wait()
            mm_tee_rc = mm_tee_proc.wait()
            mm_rc = mm_proc.wait()
            for rc, name in zip([bwa_rc, filter_rc, mm_tee_rc, mm_rc],
                                ["BWA MEM", "filter", "tee", "minimap2"]):
                if rc != 0:
                    raise typer.Exit(code=rc)


        else:
            # minimap2 subprocess
            mm_proc = subprocess.Popen(minimap2_args_list, stdout=subprocess.PIPE, stderr=log_fh)

            # filter_soft_clip subprocess
            filter_args = [
                "python", "-m", "mmbwa.filter_soft_clip.filter_script", "-",
                "--output", str(filtered_bam),
                "--output-temp", "/dev/null",
                "--threshold", str(threshold),
                "--regions_bed", str(regions_bed)]

            if unmapped:
                filter_args.append("--unmapped")

            filter_proc = subprocess.Popen(filter_args, stdin=mm_proc.stdout, stdout=subprocess.PIPE,
                                           stderr=log_fh)
            mm_proc.stdout.close()

            # bwa mem subprocess
            bwa_proc = subprocess.Popen(bwa_args_list, stdin=filter_proc.stdout, stdout=subprocess.PIPE)

            with open(realigned_bam, "wb") as bam_out:
                subprocess.run(["samtools", "view", "-b"], stdin=bwa_proc.stdout, stdout=bam_out, stderr=log_fh)

            # ensure upstream processes are reaped and check exit codes
            bwa_rc = bwa_proc.wait()
            filter_rc = filter_proc.wait()
            mm_rc = mm_proc.wait()
            for rc, name in zip([bwa_rc, filter_rc, mm_rc], ["BWA MEM", "filter", "minimap2"]):
                if rc != 0:
                    raise typer.Exit(code=rc)


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
        if not no_merge:
            if filtered_bam and filtered_bam.exists():
                filtered_bam.unlink()
                typer.echo(f"Deleted temporary file: {filtered_bam}")
        for f in (softclipped_bam, mm_sam):
            if f and f.exists():
                f.unlink()
                typer.echo(f"Deleted temporary file: {f}")

    if log_fh:
        log_fh.close()



def main():
    app()


if __name__ == "__main__":
    main()
