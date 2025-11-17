import os
import logging
from datetime import datetime
import glob
import argparse
import subprocess
import shutil


def create_log(status, random_state):
    """
    Create log file in logs directory.
    :param status: if True, create log file with current date and time; else create test.log
    :param random_state: random state used for sampling
    :return: logger object
    """
    if not os.path.exists(os.path.join(os.getcwd(), "logs")):
        os.mkdir("logs")
    if status:
        logger = logging.getLogger(__name__)
        logging.basicConfig(filename=os.path.join(os.getcwd(),
                                                  "logs",
                                                  datetime.now().strftime(
                                                      f"%Y-%m-%d_%H:%M_{__file__.split('/')[-1]}_random_state"
                                                      f"_{random_state}.log")),
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S',
                            level=logging.INFO)
    else:
        logger = logging.getLogger(__name__)
        logging.basicConfig(filename=os.path.join(os.getcwd(),
                                                  "logs",
                                                  datetime.now().strftime("test.log")),
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


def downsample(input_files, out_dir, random_state):
    """
    Downsample fastq files.
    :param input_files: list of input files
    :param out_dir: output directory
    :param random_state: random state used for sampling
    :return: downsampled files
    """
    logger.info("Starting downsample...")

    # Set downsample proportion
    downsample_proportion = {
        "S1": 1.0000, "S2": 0.2778, "S3": 0.0803, "S4": 1.0000, "S5": 0.3167,
        "S6": 0.3587, "S7": 0.0833, "S8": 1.0000, "S9": 0.8433, "S10": 0.6630,
        "S11": 0.1329, "S12": 0.1620, "S13": 1.0000, "S14": 0.2108, "S15": 0.0714,
        "S16": 1.0000, "S17": 0.3685, "S18": 0.3332, "S19": 0.0590, "S20": 1.0000,
        "S21": 0.7184, "S22": 0.6456, "S23": 0.1192, "S24": 0.1623, "S25": 1.0000
    }
    logger.info(f"Downsample proportion: {downsample_proportion}")

    # Downsample files
    for file in input_files:
        sample = file.split("/")[-1].split("_")[1]
        proportion = downsample_proportion[sample]
        if proportion == 1.0000:
            logger.info(f"Skipping {file} with proportion {proportion}...")
            continue
        logger.info(f"Downsampling {file} with proportion {proportion}...")

        # Downsample using SeqKit
        filename = file.split("/")[-1]
        out_name = f"{filename.split('_')[0]}_{filename.split('_')[1]}.fastq.gz"
        out_file = os.path.join(out_dir, out_name)
        cmd = f"seqkit sample -p {proportion} -s {random_state} {file} -o {out_file}"

        try:
            logger.info(subprocess.check_output(cmd, shell=True).decode("utf-8"))
        except subprocess.CalledProcessError as e:
            logger.error(f"Error: {e}")

    logger.info("Downsample done!")

    return None


def create_downsample_none(files, downsample_none_dir):
    downsample_proportion = {
        "S1": 1.0000, "S2": 0.2778, "S3": 0.0803, "S4": 1.0000, "S5": 0.3167,
        "S6": 0.3587, "S7": 0.0833, "S8": 1.0000, "S9": 0.8433, "S10": 0.6630,
        "S11": 0.1329, "S12": 0.1620, "S13": 1.0000, "S14": 0.2108, "S15": 0.0714,
        "S16": 1.0000, "S17": 0.3685, "S18": 0.3332, "S19": 0.0590, "S20": 1.0000,
        "S21": 0.7184, "S22": 0.6456, "S23": 0.1192, "S24": 0.1623, "S25": 1.0000
    }

    for file in files:
        sample = file.split("/")[-1].split("_")[1]
        filename = file.split("/")[-1]
        out_name = f"{filename.split('_')[0]}_{filename.split('_')[1]}.fastq.gz"
        proportion = downsample_proportion[sample]
        if proportion == 1.0000:
            # copy file to downsample_none_dir and rename it as downsample_none_dir/sample
            output_path = os.path.join(downsample_none_dir, out_name)
            shutil.copy(file, output_path)

    return None


if __name__ == "__main__":
    # Specify random state
    # Ask user to specify random state using argparse with the flag --random_state
    parser = argparse.ArgumentParser()
    parser.add_argument("--random_state", type=int, help="Random state used for sampling in SeqKit")
    args = parser.parse_args()
    random_state = args.random_state

    # Create log file
    logger = create_log(True, random_state)

    # Get root directory
    root_dir = os.getcwd()
    logger.info(f"Root directory: {root_dir}")

    # Set input directory
    in_dir = os.path.join(root_dir, "01_06_merged_and_filtered_rc")
    logger.info(f"Input directory: {in_dir}")

    # Get input files
    input_files = order_files(get_files(os.path.join(in_dir, "*filtered_rc.fastq.gz")))

    # Set output directory
    out_root_dir = create_dir(root_dir, "02_01_downsample_outputs")
    out_dir = create_dir(out_root_dir, f"downsample_{random_state}")
    logger.info(f"Output directory: {out_dir}")

    # Check if downsample_none exists in out_root_dir
    downsample_none_dir = os.path.join(out_root_dir, "downsample_none")
    if not os.path.isdir(downsample_none_dir):
        os.mkdir(downsample_none_dir)
        create_downsample_none(input_files, downsample_none_dir)


    # Downsample files
    downsample(input_files, out_dir, random_state)

    logger.info("Done!")








