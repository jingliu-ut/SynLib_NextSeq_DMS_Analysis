import os
import logging
from datetime import datetime
import glob
from Bio.Seq import Seq
import subprocess


def create_log(status):
    """
    Create log file in logs directory.
    :param status: if True, create log file with current date and time; else create test.log
    :return: logger object
    """
    if not os.path.exists(os.path.join(os.getcwd(), "logs")):
        os.mkdir("logs")

    if status:
        logger = logging.getLogger(__name__)
        log_filename = os.path.join(os.getcwd(), "logs", datetime.now().strftime(f"%Y-%m-%d_%H:%M_{__file__.split('/')[-1]}.log"))
        logging.basicConfig(filename=log_filename,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S',
                            level=logging.INFO)
    else:
        logger = logging.getLogger(__name__)
        log_filename = os.path.join(os.getcwd(), "logs", datetime.now().strftime("test.log"))
        logging.basicConfig(filename=log_filename,
                            format='%(asctime)s - [%(filename)s] - %(levelname)s - %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S',
                            level=logging.INFO)

    return logger


def create_dir(parent_dir, dir_name):
    """
    Create directory.
    :param parent_dir: parent directory
    :param dir_name: new directory name
    :return: new directory path
    """
    out_dir = os.path.join(parent_dir, dir_name)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        logger.info(f"Created directory: {out_dir}")

    return str(os.path.join(parent_dir, dir_name))


def get_files(pattern):
    """
    Get files in directory.
    :param pattern: pattern to search for files
    :return: a list of files
    """
    logger.info(f"Getting files...")
    files = sorted([file for file in glob.glob(pattern)])
    for file in files:
        logger.info(file)
    return files


def order_files(files):
    """
    Order files from 1 to n.
    :param files: files to order
    :return: ordered files
    """
    order_pattern = [f"{i+1}_" for i in range(len(files))]
    return [file for pattern in order_pattern for file in files if pattern in file
            and file.split("/")[-1].startswith(pattern)]


def get_raw_reads_files():
    """
    Get raw reads files.
    :return: forward and reverse raw reads files
    """
    # Get forward and reverse raw reads files in NextSeq2000_raw_reads directory.
    R1_raw_reads = order_files(get_files(os.path.join(raw_reads_dir, "*R1*.fastq.gz")))
    R2_raw_reads = order_files(get_files(os.path.join(raw_reads_dir, "*R2*.fastq.gz")))

    return R1_raw_reads, R2_raw_reads


def reverse_complement(primer):
    """
    Get reverse complement of primer.
    :param primer: primer sequence
    :return: reverse complement of primer
    """
    return str(Seq(primer).reverse_complement())


def cutadapt():
    """
    Cut adapters from raw reads.
    """
    logger.info("Start cutadapt...")

    # Log cutadapt path and version
    cutadapt_exe = '/home/steven/anaconda3/envs/cutnmerge/bin/cutadapt'
    logger.info(f"Cutadapt path: {subprocess.check_output(['which', 'cutadapt']).decode('utf-8').strip()}")
    logger.info(f"Cutadapt version: {subprocess.check_output([cutadapt_exe, '--version']).decode('utf-8').strip()}")

    # Create output directory.
    out_dir = create_dir(root_dir, "01_01_cutadapt")
    logger.info(f"Output directory: {out_dir}")

    # Create output file names.

    R1_cutadapt_files = [os.path.join(out_dir,
                                      f.split("/")[-1].replace(".fastq.gz", "_cutadapt.fastq.gz"))
                         for f in R1_raw_reads_files]
    logger.info(f"Forward output files:{R1_cutadapt_files}")


    R2_cutadapt_files = [os.path.join(out_dir,
                                      f.split("/")[-1].replace(".fastq.gz", "_cutadapt.fastq.gz"))
                         for f in R2_raw_reads_files]
    logger.info(f"Reverse output files:{R2_cutadapt_files}")

    # Setup cutadapt flags

    # Old paremeters with R2 as forward and R1 as reverse
    # forward_flags = ["-g", f"{R2_primer};max_errors=0.3;min_overlap=13", "-a", f"{R1_primer_rc};max_errors=0.1;min_overlap=18"]
    # reverse_flags = ["-G", f"{R1_primer};max_errors=0.2;min_overlap=15", "-A", f"{R2_primer_rc};max_errors=0.1;min_overlap=17"]

    R1_flags = ["-g", f"{R1_primer};max_errors=0.2;min_overlap=15", "-a", f"{R2_primer_rc};max_errors=0.1;min_overlap=17"]
    R2_flags = ["-G", f"{R2_primer};max_errors=0.3;min_overlap=13", "-A", f"{R1_primer_rc};max_errors=0.1;min_overlap=18"]

    cutadapt_commands = [[cutadapt_exe] +
                         R1_flags +
                         R2_flags +
                         ["-n", "2",
                         "-j", "12",
                         "--pair-filter=any",
                         "-m", "250",
                         "-M", "290",
                         "-o", R1_cutadapt_files[i],
                         "-p", R2_cutadapt_files[i],
                          R1_raw_reads_files[i],
                          R2_raw_reads_files[i]] for i in range(len(R1_raw_reads_files))]

    # Cut adapters from forward and reverse reads.
    for command in cutadapt_commands:
        logger.info(f"Running cutadapt command: {' '.join(command)}")
        try:
            logger.info(subprocess.check_output(command).decode('utf-8'))
        except subprocess.CalledProcessError as e:
            logger.error(f"Error: {e}")

    logger.info("Cutadapt done!")
    return None


def NGmerge():
    """
    Merge paired-end reads using NGmerge.
    """
    logger.info("Start NGmerge...")

    # Set output directory.
    out_dir = create_dir(root_dir, "01_02_NGmerge")
    logger.info(f"Output directory: {out_dir}")

    # Get input files.
    forward_files = order_files(get_files(os.path.join(root_dir, "01_01_cutadapt", "*R1*.fastq.gz")))
    logger.info(f"Forward files: {forward_files}")

    reverse_files = order_files(get_files(os.path.join(root_dir, "01_01_cutadapt", "*R2*.fastq.gz")))
    logger.info(f"Reverse files: {reverse_files}")

    # Create output file names.
    merged_files = [os.path.join(out_dir, f.split("/")[-1].replace("R1_001_cutadapt.fastq.gz", "merged.fastq.gz"))
                    for f in forward_files]
    logger.info(f"Output files: {merged_files}")

    # NGmerge commands.
    NGmerge_commands = [
        ["NGmerge"] + [
            "-1", forward_files[i],
            "-2", reverse_files[i],
            "-o", merged_files[i],
            "-m", "20",  # Minimum overlap of the paired-end reads
            "-p", "0.10",  # Mismatches to allow in the overlapped region (a fraction of the overlap length)
            "-s",  # Option to produce shortest stitched read
            "-j", os.path.join(out_dir, f"NGmerge_{i+1}_alignments.log"),  # Log file for formatted alignments of merged reads
            "-n", "12",  # Number of threads to use
            "-v"  # Option to print status updates/counts to stderr
        ] for i in range(len(forward_files))
    ]

    # Run NGmerge commands.
    for command in NGmerge_commands:
        logger.info(f"Running NGmerge command: {' '.join(command)}")
        try:
            logger.info(subprocess.check_output(command).decode('utf-8'))
        except subprocess.CalledProcessError as e:
            logger.error(f"Error: {e}")

    logger.info("NGmerge done!")
    return None


def length_filter(length):
    """
    Filter reads based on read length.
    :param length: minimum read length
    """
    logger.info("Start filtering reads based on read length...")

    # Set output directory.
    out_dir = create_dir(root_dir, f"01_03_filtered_reads_length_{length}")
    logger.info(f"Output directory: {out_dir}")

    # Get input files.
    merged_files = order_files(get_files(os.path.join(root_dir, "01_02_NGmerge", "*merged.fastq.gz")))
    logger.info(f"Input files: {merged_files}")

    # Create output file names.
    filtered_files = [os.path.join(out_dir, f.split("/")[-1].replace("merged.fastq.gz", "filtered.fastq.gz"))
                      for f in merged_files]
    logger.info(f"Output files: {filtered_files}")

    # Filter reads based on read length.
    filter_commands = [
        ["seqkit", "seq"] + [
            "-m", str(length),
            "-M", str(length),
            merged_files[i],
            "-o", filtered_files[i]
    ] for i in range(len(merged_files))]

    # Run filter commands.
    for command in filter_commands:
        logger.info(f"Running filter command: {' '.join(command)}")
        try:
            logger.info(subprocess.check_output(command).decode('utf-8'))
        except subprocess.CalledProcessError as e:
            logger.error(f"Error: {e}")

    logger.info("Filtering done!")
    return None


def qc_filter(qc):
    """
    Filter reads based on read quality.
    :param length: minimum read quality
    """
    logger.info("Start filtering reads based on read quality...")

    # Set output directory.
    out_dir = create_dir(root_dir, f"01_04_filtered_reads_qc_{qc}")
    logger.info(f"Output directory: {out_dir}")

    # Get input files.
    length_files = order_files(get_files(os.path.join(root_dir, "01_03_filtered_reads_length_303", "*filtered.fastq.gz")))
    logger.info(f"Input files: {length_files}")

    # Create output file names.
    qc_files = [os.path.join(out_dir, f.split("/")[-1])
                      for f in length_files]
    logger.info(f"Output files: {qc_files}")

    # Filter reads based on read length.
    filter_commands = [
        ["fastp"] + [
            "-A",
            "-e", str(qc),
            "-i", length_files[i],
            "-o", qc_files[i]
    ] for i in range(len(length_files))]

    # Run filter commands.
    for command in filter_commands:
        logger.info(f"Running filter command: {' '.join(command)}")
        try:
            logger.info(subprocess.check_output(command).decode('utf-8'))
        except subprocess.CalledProcessError as e:
            logger.error(f"Error: {e}")

    logger.info("Filtering done!")
    return None


def removeN():
    """
    Remove reads with ambiguous bases.
    :return: None
    """
    logger.info("Start filtering reads with ambiguous bases...")

    # Set output directory.
    out_dir = create_dir(root_dir, f"01_05_filtered_reads_N")
    logger.info(f"Output directory: {out_dir}")

    # Get input files.
    qc_files = order_files(get_files(os.path.join(root_dir, "01_04_filtered_reads_qc_33", "*filtered.fastq.gz")))
    logger.info(f"Input files: {qc_files}")

    # Create output file names.
    n_files = [os.path.join(out_dir, f.split("/")[-1])
                for f in qc_files]
    logger.info(f"Output files: {n_files}")

    # Filter reads based on read length.
    filter_commands = [
        ["fastp"] + [
            "-A",
            "-n", str(0),
            "-i", qc_files[i],
            "-o", n_files[i]
        ] for i in range(len(qc_files))]

    # Run filter commands.
    for command in filter_commands:
        logger.info(f"Running filter command: {' '.join(command)}")
        try:
            logger.info(subprocess.check_output(command).decode('utf-8'))
        except subprocess.CalledProcessError as e:
            logger.error(f"Error: {e}")

    logger.info("Filtering done!")
    return None


def reverse_complement_reads():
    """
    Reverse complement reads.
    :return: None
    """
    logger.info("Start reverse complementing reads...")

    # Set output directory.
    out_dir = create_dir(root_dir, f"01_06_merged_and_filtered_rc")
    logger.info(f"Output directory: {out_dir}")

    # Get input files.
    in_files = order_files(get_files(os.path.join(root_dir, "01_05_filtered_reads_N", "*filtered.fastq.gz")))
    logger.info(f"Input files: {in_files}")

    # Create output file names.
    rc_files = [os.path.join(out_dir, f.split("/")[-1].replace("filtered.fastq.gz", "filtered_rc.fastq.gz")) for f in in_files]
    logger.info(f"Output files: {rc_files}")

    # Reverse complement reads.
    rc_commands = [
        ["seqkit", "seq"] + [
            "-r", "-p",
            "-t", "DNA",
            in_files[i],
            "-o", rc_files[i]
        ] for i in range(len(in_files))]

    # Run reverse complement commands.
    for command in rc_commands:
        logger.info(f"Running reverse complement command: {' '.join(command)}")
        try:
            logger.info(subprocess.check_output(command).decode('utf-8'))
        except subprocess.CalledProcessError as e:
            logger.error(f"Error: {e}")

    logger.info("Reverse complement done!")
    return None


if __name__ == "__main__":
    # Create log file.
    logger = create_log(True)

    # Get root directory.
    root_dir = os.getcwd()
    logger.info(f"Root directory: {root_dir}")
    raw_reads_dir = os.path.join(root_dir, "00_NextSeq2000_raw_reads")

    # Get forward and reverse reads files and store in separate lists.
    R1_raw_reads_files, R2_raw_reads_files = get_raw_reads_files()
    R1_raw_reads_files = order_files(R1_raw_reads_files)
    R2_raw_reads_files = order_files(R2_raw_reads_files)

    # Set forward and reverse primers and get their reverse complements.
    R2_primer = 'GCCGCACCCCCACTCGCCGGCTGGTCC'
    R1_primer = 'CCGAAGTTGGAGCCCTGGTG'
    R2_primer_rc = reverse_complement(R2_primer)
    R1_primer_rc = reverse_complement(R1_primer)
    logger.info(f"Forward primer: {R2_primer}")
    logger.info(f"Reverse primer: {R1_primer}")
    logger.info(f"Forward primer reverse complement: {R2_primer_rc}")
    logger.info(f"Reverse primer reverse complement: {R1_primer_rc}")

    # Run cutadapt.
    cutadapt()

    # Run NGmerge.
    NGmerge()

    # Filter reads based on read length.
    length_filter(length=303)

    # Filter reads based on quality score.
    qc_filter(qc=33)

    # Ambiguous base filtering
    removeN()

    # Reverse complement reads
    reverse_complement_reads()

    logger.info("Done!")