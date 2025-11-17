import glob
from datetime import datetime
import os
import logging
import vcfpy
import gzip
import argparse


def create_log(status):
    """
    Create log file in logs directory.
    :param status: if True, create log file with current date and time; else create test.log
    :return: logger object
    """
    if not os.path.exists(os.path.join(os.getcwd(), "logs")):
        os.mkdir("logs")

    sample = get_sample_name(input_file)

    if status:
        logger = logging.getLogger(__name__)
        log_filename = os.path.join(os.getcwd(), "logs", datetime.now().strftime(f"%Y-%m-%d_%H:%M_{__file__.split('/')[-1]}_downsample{downsample_num}_{sample}.log"))
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


def process_string(seq):
    """
    Read fastq file.
    :param file: fastq file
    :return: list of reads and quality scores
    """
    seq_codon = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    seq_dict = {}

    for i in range(len(seq_codon)):
        seq_dict[i + 177] = seq_codon[i]

    return seq_dict


def get_sample_name(file):
    """
    Get sample name from id.
    :param id: id
    :return: sample name
    """
    sample_dict = {
        "S1": "R1_G1", "S2": "R1_G2", "S3": "R1_G3", "S4": "R1_L1", "S5": "R1_L2",
        "S6": "R1_L3", "S7": "R1_L4", "S8": "R1_Y1", "S9": "R1_Y2", "S10": "R1_Y3",
        "S11": "R1_Y4", "S12": "R1_Y5", "S13": "R2_G1", "S14": "R2_G2", "S15": "R2_G3",
        "S16": "R2_L1", "S17": "R2_L2", "S18": "R2_L3", "S19": "R2_L4", "S20": "R2_Y1",
        "S21": "R2_Y2", "S22": "R2_Y3", "S23": "R2_Y4", "S24": "R2_Y5", "S25": "WT"
    }
    id = file.split("/")[-1].split("_")[1].split(".")[0]

    return sample_dict[id]


def check_variant(seq_dict, qual_dict, ref_dict_nt, ref_dict_aa, writer, variant_id, read_id, replicate, sample_name):
    """
    Check for variants and write to vcf file.
    :param seq_dict: sequence in dictionary format
    :param qual_dict: quality scores in dictionary format
    :param ref_dict_nt: reference nucleotide in dictionary format
    :param ref_dict_aa: reference amino acid in dictionary format
    :param writer: writer object
    :param variant_id: variant id
    :param read_id: read id
    :param replicate: replicate number
    :param sample_name: sample name
    :return: variant count
    """

    variant_count = 0

    for codon in seq_dict:
        ref_nt = ref_dict_nt[codon]
        alt_nt = seq_dict[codon]
        if ref_nt != alt_nt:
            variant_count += 1
            data = {
                "variant_id": variant_id + "_" + str(variant_count),
                "read_id": read_id,
                "replicate": replicate,
                "sample_name": sample_name,
                "POS_NT": codon*3-2,
                "REF_NT": ref_nt,
                "ALT_NT": alt_nt,
                "REF_AA": ref_dict_aa[codon],
                "ALT_AA": translate_dna(alt_nt),
                "AA_POS": codon,
                "AA_CHANGE": ref_dict_aa[codon] + str(codon) + translate_dna(alt_nt),
                "TYPE": check_type(ref_nt, alt_nt),
                "BQ": gen_BQ(qual_dict[codon])
            }
            vcf_record(writer, data)



    return variant_count


def variant_calling(file, output_dir):
    """
    Call variants from fastq file and write to vcf file.
    :param file: fastq file
    :param output_dir: output directory
    :return: None
    """

    logger.info(f"Calling variants from {file}...")

    ref_dict_nt, ref_dict_aa = gen_ref_dict()

    sample = get_sample_name(file)
    out_file = os.path.join(output_dir, f"{sample}.vcf")

    # Create vcf file
    writer = create_vcf(out_file)
    logger.info(f"Writing to {out_file}...")


    seq_count = 0
    qual_count = 0
    variant_count = 0


    file = gzip.open(file, "rt")
    with file as f:
        for line in f:
            if line.startswith("@"):
                seq = f.readline().strip()
                seq_dict = process_string(seq)
                seq_count += 1
                f.readline()
                qual = f.readline().strip()
                qual_dict = process_string(qual)
                qual_count += 1

                variant_id = sample + "_" + str(seq_count)
                read_id = seq_count
                sample_name = None
                if sample == "WT":
                    sample_name = "WT"
                    replicate = 0
                else:
                    sample_name = sample.split("_")[1]
                    replicate = int(sample[1])
                variant_count += check_variant(seq_dict, qual_dict, ref_dict_nt, ref_dict_aa, writer,
                                               variant_id, read_id, replicate, sample_name)

    logger.info(f"Number of reads loaded: {seq_count}")
    if seq_count == qual_count:
        logger.info("Number of reads and quality scores match.")

    logger.info(f"Total number of variants called: {variant_count}")

    return None


def gen_BQ(qual):
    """
    Generate base quality.
    :param qual: quality scores
    :return: base quality
    """
    bq = []
    for i in qual:
        bq.append(qc_converter(i))

    return bq


def check_type(ref, alt):
    """
    Check if variant is SNP or MNP.
    :param ref: reference nucleotide
    :param alt: alternate nucleotide
    :return: SNP or MNP
    """
    mismatch = 0
    for i in range(len(ref)):
        if ref[i] != alt[i]:
            mismatch += 1

    if mismatch == 1:
        return "SNP"
    else:
        return "MNP"


def gen_ref_dict():
    seq = "AGGTACATCCCCGAGGGCCTGCAGTGCTCGTGTGGAATCGACTACTACACGCTCAAGCCGGAGGTCAACAACGAGTCTTTTGTCATCTACATGTTCGTGGTCCACTTCACCATCCCCATGATTATCATCTTTTTCTGCTATGGGCAGCTCGTCTTCACCGTCAAGGAGGCCGCTGCCCAGCAGCAGGAGTCAGCCACCACACAGAAGGCAGAGAAGGAGGTCACCCGCATGGTCATCATCATGGTCATCGCTTTCCTGATCTGCTGGGTGCCCTACGCCAGCGTGGCATTCTACATCTTCACC"
    seq_codon = [seq[i:i+3] for i in range(0, len(seq), 3)]
    ref_dict_nt = {}

    for i in range(len(seq_codon)):
        ref_dict_nt[i+177] = seq_codon[i]

    ref_dict_aa = {}
    for i in ref_dict_nt:
        ref_dict_aa[i] = translate_dna(ref_dict_nt[i])

    return ref_dict_nt, ref_dict_aa


def qc_converter(qual):
    """
    Convert quality scores to base quality.
    :param qual: quality scores
    :return: base quality
    """
    return ord(qual) - 33


def create_vcf(out_name):
    """
    Create vcf file.
    :param out_name: output file name
    :return: writer object
    """
    header = vcfpy.Header(
        lines=[
            vcfpy.HeaderLine(key='fileformat', value='VCFv4.3'),
            vcfpy.FilterHeaderLine.from_mapping(mapping={'ID': '.', 'Description': 'no filter'}),
            vcfpy.InfoHeaderLine.from_mapping(
                mapping={'ID': 'variant_id', 'Number': '1', 'Type': 'String',
                         'Description': '[sample name]_[read number]_[variant number]'}),
            vcfpy.InfoHeaderLine.from_mapping(
                mapping={'ID': 'read_id', 'Number': '1', 'Type': 'Integer',
                         'Description': 'read number'}
            ),
            vcfpy.InfoHeaderLine.from_mapping(
                mapping={'ID': 'replicate', 'Number': '1', 'Type': 'Integer',
                         'Description': 'replicate number'}
            ),
            vcfpy.InfoHeaderLine.from_mapping(
                mapping={'ID': 'sample', 'Number': '1', 'Type': 'String',
                         'Description': 'sample name'}
            ),
            vcfpy.InfoHeaderLine.from_mapping(
                mapping={'ID': 'POS_NT', 'Number': '1', 'Type': 'Integer', 'Description': 'Nucleotide position'}),
            vcfpy.InfoHeaderLine.from_mapping(
                mapping={'ID': 'REF_NT', 'Number': '1', 'Type': 'String', 'Description': 'Reference nucleotide'}),
            vcfpy.InfoHeaderLine.from_mapping(
                mapping={'ID': 'ALT_NT', 'Number': '1', 'Type': 'String', 'Description': 'Alternate nucleotide'}),
            vcfpy.InfoHeaderLine.from_mapping(
                mapping={'ID': 'REF_AA', 'Number': '1', 'Type': 'String', 'Description': 'Reference amino acid'}),
            vcfpy.InfoHeaderLine.from_mapping(
                mapping={'ID': 'ALT_AA', 'Number': '1', 'Type': 'String', 'Description': 'Alternate amino acid'}),
            vcfpy.InfoHeaderLine.from_mapping(
                mapping={'ID': 'AA_CHANGE', 'Number': '1', 'Type': 'String',
                         'Description': '[ref AA][pos AA][alt AA]'}),
            vcfpy.InfoHeaderLine.from_mapping(
                mapping={'ID': 'AA_POS', 'Number': '1', 'Type': 'Integer', 'Description': 'Amino acid position'}),
            vcfpy.InfoHeaderLine.from_mapping(
                mapping={'ID': 'TYPE', 'Number': '1', 'Type': 'String', 'Description': 'SNP or MNP'}),
            vcfpy.InfoHeaderLine.from_mapping(
                mapping={'ID': 'BQ', 'Number': '3', 'Type': 'Integer', 'Description': 'Base quality'})
        ],
        samples=vcfpy.SamplesInfos([])
    )

    writer = vcfpy.Writer.from_path(out_name, header)

    return writer


def vcf_record(writer, data):
    """
    Write record to vcf file.
    :param writer:
    :param data:
    :return:
    """
    record = vcfpy.Record(
        CHROM='HuRHWT_ymS',
        POS=data['POS_NT'],
        ID=".",
        REF=data['REF_NT'],
        ALT=[vcfpy.Substitution(type_=data['TYPE'], value=data['ALT_NT'])],
        QUAL=".",
        FILTER=".",
        INFO={"variant_id": data['variant_id'],
              "read_id": data['read_id'],
              "replicate": data['replicate'],
              "sample": data['sample_name'],
              "POS_NT": data['POS_NT'],
              "REF_NT": data['REF_NT'],
              "ALT_NT": data['ALT_NT'],
              "REF_AA": data['REF_AA'],
              "ALT_AA": data['ALT_AA'],
              "AA_POS": data['AA_POS'],
              "AA_CHANGE": data['AA_CHANGE'],
              "TYPE": data['TYPE'],
              "BQ": data['BQ']
              },
        FORMAT=[],
        calls=[]
    )

    writer.write_record(record)

    return None


def translate_dna(dna):
    """
    Translate single DNA sequence to amino acid sequence.
    :param dna: dna sequence
    :return: amino acid sequence
    """
    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }
    aa = ""
    for i in range(0, len(dna), 3):
        codon = dna[i:i + 3]
        if len(codon) == 3:
            aa += codon_table[codon]
    return aa


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--random_state", type=int, help="Random state used for sampling in SeqKit")
    parser.add_argument("-i", "--input_file", type=str, help="Input file")
    args = parser.parse_args()
    downsample_num = args.random_state
    input_file = args.input_file

    if downsample_num == -1:
        downsample_num = "none"

    # Create log file.
    logger = create_log(True)

    # Get root directory.
    root_dir = os.getcwd()
    logger.info(f"Root directory: {root_dir}")


    # Set output directory.
    downsample_output_root = create_dir(root_dir, "03_01_vcf_files")
    logger.info(f"Output root directory: {downsample_output_root}")
    downsample_output_dir = create_dir(downsample_output_root, f"downsample_{downsample_num}")
    logger.info(f"Output directory: {downsample_output_dir}")

    # Get files in input directory.
    logger.info(f"Input file: {input_file}")

    try:
        variant_calling(input_file, downsample_output_dir)
    except Exception as e:
        logger.error(f"Error: {e}")


    logger.info("Done!")