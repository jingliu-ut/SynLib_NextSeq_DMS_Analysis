import allel
import glob
import pandas as pd
from datetime import datetime
import argparse
import os
import logging
import cudf
import matplotlib.pyplot as plt
import numpy as np
import dask.dataframe as dd


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


def get_undownsample_files(pattern, root_dir):
    """
    Get files in directory for undownsample samples.
    :param pattern: pattern to search for files
    :return: a list of files
    """
    logger.info(f"Getting files...")
    undownsample_samples = [
        "R1_G1", "R1_L1", "R1_Y1",
        "R2_G1", "R2_L1", "R2_Y1",
        "WT"
    ]
    undownsample_dir = os.path.join(root_dir, "downsample_none")
    pattern = os.path.join(undownsample_dir, pattern)
    files = sorted([file for file in glob.glob(pattern) if file.split("/")[-1].split(".")[0] in undownsample_samples])
    return files


def vcf_to_csv(file, out_dir):
    """
    Convert VCF file to CSV file.
    :param file: VCF file
    :return: None
    """
    sample = file.split("/")[-1].split(".")[0]

    out_file = os.path.join(out_dir, f"{sample}.csv")

    logger.info(f"Converting VCF to CSV: {file}")
    allel.vcf_to_csv(
        file,
        out_file,
        fields=["read_id", "replicate", "sample", "POS_NT", "REF_NT", "ALT_NT", "AA_CHANGE", "BQ"]
    )

    return None


def concat_replicates(files):
    bins = ["G1", "G2", "G3", "L1", "L2", "L3", "L4", "Y1", "Y2", "Y3", "Y4", "Y5"]

    for bin in bins:
        temp_files = [file for file in files if bin in file]

        # Concatenate R1 and R2 csv files
        ddf = dd.concat([dd.read_csv(file) for file in temp_files], ignore_index=True)

        # Save concatenated csv to downsample directory
        out_dir = create_dir(out_downsample_dir, "raw_csv_concat")
        out_file = os.path.join(out_dir, f"{bin}.parquet")
        ddf.to_parquet(out_file)


    return None



def bq_filter(file, bq_threshold, out_dir):
    """
    Filter out variants with BQ < 30.
    :param file: VCF file
    :return: None
    """
    logger.info(f"Applying BQ filter with BQ >= {bq_threshold}...")

    sample = file.split("/")[-1].split(".")[0]

    # Read csv into cudf
    raw_ddf = dd.read_csv(file)
    pre_filter_len = len(raw_ddf)
    logger.info(f"Number of reads before bq filter: {pre_filter_len}")

    # Filter out variants with BQ_1 < bq_threshold or BQ_2 < bq_threshold or BQ_3 < bq_threshold
    filtered_ddf = raw_ddf[(raw_ddf["BQ_1"] >= bq_threshold) & (raw_ddf["BQ_2"] >= bq_threshold) & (raw_ddf["BQ_3"] >= bq_threshold)]
    post_filter_len = len(filtered_ddf)
    logger.info(f"Number of reads after bq filter: {post_filter_len}")
    logger.info(f"Percentage of reads kept: {post_filter_len/pre_filter_len*100:.2f}%")

    # Save variant_id to csv
    out_file = os.path.join(out_dir, f"{sample}.parquet")
    filtered_ddf.to_parquet(out_file, engine="pyarrow")

    return None


def get_read_id(df):
    """
    Extract read_id from variant_id.
    :param df: dataframe
    :return: dataframe with read_id column
    """
    # Extract replicate number from variant_id. Split variant_id by _ and get first element, then get the second element from the string.
    df["replicate"] = df["variant_id"].str.split("_", expand=True)[0].str[1].astype("int")
    df["read_id"] = df["variant_id"].str.split("_", expand=True)[2].astype("int")
    return df


def check_type(df, sample):
    """
    Generate summaries for each sample.
    :param file:
    :return: None
    """

    # SNP and MNP summary
    # Read csv into cudf with only variant_id and TYPE columns
    type_df = df[["variant_id", "TYPE"]]

    # Count number of SNPs and MNPs for sample
    snp_count = len(type_df[type_df["TYPE"] == "SNP"])
    mnp_count = len(type_df[type_df["TYPE"] == "MNP"])
    num_reads = snp_count + mnp_count

    # Extract read_id from variant_id
    get_read_id(type_df)

    # Organize date by read_id and TYPE and count number of each TYPE for each read_id. Pivot TYPE column to SNP and
    # MNP columns and fill with count.
    type_df = type_df.groupby(["read_id", "TYPE"]).size().reset_index(name="count")
    type_df = type_df.pivot(index="read_id", columns="TYPE", values="count").reset_index().fillna(0)
    # Reorder columns
    type_df = type_df[["read_id", "SNP", "MNP"]]

    # Add snp_freq and mnp_freq columns to read_df for each read_id. snp_freq = snp_count/num_reads for each read_id.
    # Add snp_freq and mnp_freq and total columns to type_df
    type_df["total"] = type_df["SNP"] + type_df["MNP"]
    # Collapse type_df to get unique read_id and count for each read_id
    type_df["snp_freq"] = type_df["SNP"] / type_df["total"]
    type_df["mnp_freq"] = type_df["MNP"] / type_df["total"]

    total_df = type_df["total"].value_counts().reset_index().to_pandas()

    plt.scatter(total_df["total"], total_df["count"], s=8, alpha=1, c=total_df["total"], cmap="cool")
    plt.show()
    return None


def check_abundance(df, sample):



    variant_df = df[["POS_NT", "REF_NT", "ALT_NT", "AA_CHANGE", "AA_POS"]]
    variant_df["NT_CHANGE"] = variant_df["REF_NT"] + variant_df["POS_NT"].astype("str") + variant_df["ALT_NT"]
    variant_df.drop(["REF_NT", "ALT_NT"], axis=1, inplace=True)

    # Count number of unique variant (NT_CHANGE)
    variant_count_df = variant_df.groupby("NT_CHANGE").size().reset_index(name="count")
    count_distribution = variant_count_df.groupby("count").size().reset_index(name="num_variants").to_pandas()

    # Create scatter plot for count_distribution
    plt.scatter(count_distribution["count"], count_distribution["num_variants"], s=8, alpha=1, c=count_distribution["num_variants"], cmap="cool")
    plt.xlim(-10, 100)
    plt.title(f"Abundance Distribution ({sample})")
    plt.xlabel("Abundance")
    plt.ylabel("Number of Unique Variants")
    plt.show()
    return None


def check_position(df, sample):



    variant_df = df[["variant_id", "AA_POS", "TYPE"]]

    # Count number of variants at each position and create a new column "num_variants" for each position
    position_df = variant_df.groupby(["AA_POS", "TYPE"]).size().reset_index(name="num_variants").to_pandas()

    # Create scatter plot for position_df, color by TYPE    plt.scatter(position_df["AA_POS"], position_df["num_variants"]
    # plt SNP. Set colour to #25C3DA
    plt.scatter(x="AA_POS", y="num_variants", data=position_df[position_df["TYPE"] == "SNP"], s=8, alpha=1, c="royalblue", label="SNP")
    # plt MNP. Set colour to #FF6347
    plt.scatter(x="AA_POS", y="num_variants", data=position_df[position_df["TYPE"] == "MNP"], s=8, alpha=1, c="tomato", label="MNP")

    # plt.ylim(-10, 100000)

    plt.title(f"Variant Position Distribution ({sample})")
    plt.xlabel("AA Position")
    plt.ylabel("Number of Variants")
    plt.legend(["SNP", "MNP"], loc="upper right", frameon=False)
    plt.show()

    return None


def get_single_double_variants(file, out_dir):

    logger.info(f"Getting single and double variants for {file}...")

    sample = file.split("/")[-1].split(".")[0]

    variant_ddf = dd.read_parquet(file)

    variant_ddf = variant_ddf[["read_id", "replicate", "REF_NT", "ALT_NT", "POS_NT"]]

    def combine_columns(df):
        df["NT_CHANGE"] = df["REF_NT"] + df["POS_NT"].astype(str) + df["ALT_NT"]
        return df

    variant_ddf = variant_ddf.map_partitions(combine_columns)

    variant_ddf = variant_ddf[["read_id", "replicate", "NT_CHANGE"]]

    per_read_variant_count = variant_ddf.groupby(["read_id", "replicate"]).count().reset_index()
    single_double_id = per_read_variant_count[per_read_variant_count["NT_CHANGE"] <= 2][["read_id", "replicate"]]
    del per_read_variant_count

    variant_ddf = variant_ddf.merge(single_double_id, on=["read_id", "replicate"], how="inner")

    variant_df = variant_ddf.compute()
    del variant_ddf

    variant_df["variant_num"] = variant_df.groupby(["read_id", "replicate"]).cumcount() + 1
    variant_df = variant_df.pivot(index=["read_id", "replicate"], columns="variant_num", values="NT_CHANGE").fillna(
        "").reset_index()
    variant_df.columns = ["read_id", "replicate", "variant_1", "variant_2"]

    variant_df["NT_CHANGE"] = variant_df["variant_1"] + "_" + variant_df["variant_2"]
    variant_df = variant_df.drop(["variant_1", "variant_2"], axis=1)

    variant_df["NT_CHANGE"] = variant_df["NT_CHANGE"].str.rstrip("_")

    variant_ddf = dd.from_pandas(variant_df)
    del variant_df

    variant_ddf = variant_ddf[["NT_CHANGE"]]
    variant_ddf = variant_ddf.groupby("NT_CHANGE").size().reset_index()
    variant_ddf = variant_ddf.rename(columns={0: "count"})

    logger.info(f"Number of unique variants: {len(variant_ddf['NT_CHANGE'].unique())}")

    logger.info(f"Writing single and double variants to {out_dir}/{sample}.parquet")
    variant_ddf.to_parquet(os.path.join(out_dir, f"{sample}.parquet"))
    logger.info(f"Finished writing single and double variants to {out_dir}/{sample}.parquet")

    return None


def translate(x):
    """
    Translate nucleotide change to amino acid change.
    :param x: Nucleotide change
    :return: Amino acid change
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
    variants = x.split("_")
    list_aa = [codon_table[variant[:3]] + str(int((int(variant[3:6])+2)/3)) + codon_table[variant[6:]] for variant in variants]
    aa = "_".join(list_aa)

    return aa


def merge_df(files, replicate):
    """
    Merge dataframes for each replicate.
    :param files:
    :param replicate:
    :return:
    """

    logger.info(f"Merging {files} ...")
    files = [file for file in files if replicate in file]
    files = sorted(files, key=lambda x: x.split("/")[-1])
    merged_df = pd.DataFrame(columns=["NT_CHANGE"])
    merged_ddf = dd.from_pandas(merged_df, npartitions=1)
    del merged_df
    for file in files:
        ddf = dd.read_parquet(file)
        sample = file.split("/")[-1].split(".")[0].split("_")[1]
        ddf = ddf.rename(columns={"count": sample})
        merged_ddf = dd.merge(merged_ddf, ddf, on="NT_CHANGE", how="outer").fillna(0)
        del ddf
    merged_ddf["total"] = merged_ddf.iloc[:, 1:].sum(axis=1)

    merged_ddf = merged_ddf.assign(num_variants=merged_ddf["NT_CHANGE"].str.count('_') + 1).query("num_variants <= 2").drop("num_variants", axis=1)
    merged_df = merged_ddf.compute()
    merged_df['AA_CHANGE'] = merged_df['NT_CHANGE'].apply(lambda x: translate(x))

    return merged_df


def filter_df(df, threshold, sample, replicate):

    logger.info(f"Applying abundance filter with threshold >= {threshold} for replicate: {replicate} sample: {sample}...")
    df = df[df["total"] >= threshold]
    df = df.drop("NT_CHANGE", axis=1)
    first_col = df.pop("AA_CHANGE")
    df.insert(0, "AA_CHANGE", first_col)

    # Combine duplicate rows and sum up the counts
    df = df.groupby("AA_CHANGE").sum().reset_index()

    out_dir = create_dir(out_downsample_dir, f"04_02_01_filtered_abundance_{threshold}")
    df.to_csv(os.path.join(out_dir, f"{replicate}_{sample}.csv"), index=False)

    return df



def binned_distribution(file, bin):
    """
    Generate binned distribution for each replicate and bin.
    :param merged_df: Merged dataframe
    :param bin: Bin name
    :param replicate: Replicate name
    :return: None
    """

    merged_df = pd.read_csv(file)

    # Calculate the frequency for each variant
    for column in merged_df.columns[1:-1]:
        merged_df[column] = merged_df[column] / merged_df["total"]
    # Drop "total" column
    merged_df.drop("total", axis=1, inplace=True)
    merged_df = merged_df.rename(columns={"AA_CHANGE": ""})

    out_dir = create_dir(root_dir, "04_02_dSortSeq_inputs")
    out_downsample_dir = create_dir(out_dir, f"downsample_{downsample_num}")

    if bin == "R1_G":
        # Export merged_df to csv with header and without index
        merged_df.to_csv(os.path.join(out_downsample_dir, f"R1_dark_binned_distribution.csv"), index=False)

    elif bin == "R2_G":
        # Export merged_df to csv with header and without index
        merged_df.to_csv(os.path.join(out_downsample_dir, f"R2_dark_binned_distribution.csv"), index=False)

    elif bin == "R1_L":
        # Export merged_df to csv with header and without index
        merged_df.to_csv(os.path.join(out_downsample_dir, f"R1_light_binned_distribution.csv"), index=False)

    elif bin == "R2_L":
        # Export merged_df to csv with header and without index
        merged_df.to_csv(os.path.join(out_downsample_dir, f"R2_light_binned_distribution.csv"), index=False)

    elif bin == "R1_Y":
        # Export merged_df to csv with header and without index
        merged_df.to_csv(os.path.join(out_downsample_dir, f"R1_ymScarlet_binned_distribution.csv"), index=False)

    elif bin == "R2_Y":
        # Export merged_df to csv with header and without index
        merged_df.to_csv(os.path.join(out_downsample_dir, f"R2_ymScarlet_binned_distribution.csv"), index=False)


    return None


def mixing_coefficient(file, bin):
    """
    Calculate mixing coefficient: the proportion of each variant in the library.
    :param merged_df: Merged dataframe
    :param bin: Bin name
    :param replicate: Replicate name
    :return: None
    """

    merged_df = pd.read_csv(file)
    # Calculate mixing coefficient
    total_variants = sum(merged_df["total"])
    merged_df["mixing_coefficient"] = merged_df["total"] / total_variants
    merged_df = merged_df.rename(columns={"AA_CHANGE": ""})
    merged_df = merged_df[["", "mixing_coefficient"]]

    out_downsample_dir = os.path.join(root_dir, "04_02_dSortSeq_inputs", f"downsample_{downsample_num}")

    if bin == "R1_G":
        # Export merged_df to csv without header and without index
        merged_df.to_csv(os.path.join(out_downsample_dir, f"R1_dark_mixing_coefficient.csv"), header=False,
                         index=False)
    elif bin == "R2_G":
        # Export merged_df to csv without header and without index
        merged_df.to_csv(os.path.join(out_downsample_dir, f"R2_dark_mixing_coefficient.csv"), header=False,
                            index=False)

    elif bin == "R1_L":
        # Export merged_df to csv without header and without index
        merged_df.to_csv(os.path.join(out_downsample_dir, f"R1_light_mixing_coefficient.csv"), header=False,
                         index=False)

    elif bin == "R2_L":
        # Export merged_df to csv without header and without index
        merged_df.to_csv(os.path.join(out_downsample_dir, f"R2_light_mixing_coefficient.csv"), header=False,
                         index=False)

    elif bin == "R1_Y":
        # Export merged_df to csv without header and without index
        merged_df.to_csv(os.path.join(out_downsample_dir, f"R1_ymScarlet_mixing_coefficient.csv"), header=False,
                         index=False)

    elif bin == "R2_Y":
        # Export merged_df to csv without header and without index
        merged_df.to_csv(os.path.join(out_downsample_dir, f"R2_ymScarlet_mixing_coefficient.csv"), header=False,
                         index=False)
    return None


def sorting_boundaries(flow_csv_files):
    """
    Generate sorting boundaries for each replicate and bin.
    :param flow_csv_files: Flow cytometry data files
    :return: None
    """
    for file in flow_csv_files:
        if "LIB_L.csv" in file:
            light = cudf.read_csv(file, usecols=["Comp-B510-A"]).rename(columns={"Comp-B510-A": "GFP"}).sort_values(by="GFP").to_pandas()
            ymScarlet = cudf.read_csv(file, usecols=["Comp-YG602-A"]).rename(columns={"Comp-YG602-A": "ymScarlet"}).sort_values(by="ymScarlet").to_pandas()

        elif "LIB_D.csv" in file:
            dark = cudf.read_csv(file, usecols=["Comp-B510-A"]).rename(columns={"Comp-B510-A": "GFP"}).sort_values(by="GFP").to_pandas()

    # Drop negative values
    light = light[light["GFP"] >= 0]
    dark = dark[dark["GFP"] >= 0]
    ymScarlet = ymScarlet[ymScarlet["ymScarlet"] >= 0]

    # Generate distribution plots for each dataframe and log transform values in the dataframe. If the value is negative, -log10(abs(x)) is used.
    log_light = light.map(lambda x: x if x ==0 else np.log10(x))
    log_dark = dark.map(lambda x: x if x ==0 else np.log10(x))
    log_ymScarlet = ymScarlet.map(lambda x: x if x ==0 else np.log10(x))

    # Remove nan values from the dataframe
    log_light = log_light.dropna()
    log_dark = log_dark.dropna()
    log_ymScarlet = log_ymScarlet.dropna()

    out_downsample_dir = os.path.join(root_dir, "04_02_dSortSeq_inputs", f"downsample_{downsample_num}")

    # Export the log transformed dataframes to csv
    log_light.to_csv(os.path.join(out_downsample_dir, "light_overall_distribution.csv"), header=False, index=False)
    log_dark.to_csv(os.path.join(out_downsample_dir, "dark_overall_distribution.csv"), header=False, index=False)
    log_ymScarlet.to_csv(os.path.join(out_downsample_dir, "ymScarlet_overall_distribution.csv"), header=False, index=False)

    gates = {
        "G1": 72.1, "G2": 94.5,
        "L1": 52.8, "L2": 73.5, "L3": 96,
        "Y1": 30.7, "Y2": 60.9, "Y3": 89.3, "Y4": 94.7
    }

    # Calculate the sorting boundaries for each condition
    light_boundary = pd.DataFrame(index=["L1", "L2", "L3"], columns=[""])
    dark_boundary = pd.DataFrame(index=["G1", "G2"], columns=[""])
    ymScarlet_boundary = pd.DataFrame(index=["Y1", "Y2", "Y3", "Y4"], columns=[""])

    # Find the value for the proportion of each gate
    for gate in gates:
        ind = int(gate[1]) - 1
        if gate[0] == "G":
            dark_boundary.iloc[ind] = log_dark.quantile(gates[gate]/100)[0]


        elif gate[0] == "L":
            light_boundary.iloc[ind] = log_light.quantile(gates[gate]/100)[0]


        elif gate[0] == "Y":
            ymScarlet_boundary.iloc[ind] = log_ymScarlet.quantile(gates[gate]/100)[0]


    # export the dark_boundary to csv
    dark_boundary.to_csv(os.path.join(out_downsample_dir, "dark_boundaries.csv"), header=False, index=True)
    # export the light_boundary to csv
    light_boundary.to_csv(os.path.join(out_downsample_dir, "light_boundaries.csv"), header=False, index=True)
    # export the ymScarlet_boundary to csv
    ymScarlet_boundary.to_csv(os.path.join(out_downsample_dir, "ymScarlet_boundaries.csv"), header=False, index=True)

    return None


if __name__ == '__main__':

    # Create log file.
    logger = create_log(True)

    # Get root directory.
    root_dir = os.getcwd()
    logger.info(f"Root directory: {root_dir}")


    # Get downsample folder number.
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--number", help="Downsample folder number", type=int, required=True)
    parser.add_argument("-bq", "--bq", help="BQ threshold", type=int, default=33)
    parser.add_argument("-a", "--abundance", help="Abundance threshold", type=int, default=50)
    args = parser.parse_args()
    bq_threshold = args.bq
    downsample_num = args.number
    threshold = args.abundance
    downsample_dir = os.path.join("04_01_variant_analysis", f"downsample_{downsample_num}")
    logger.info(f"Downsample directory: {downsample_dir}")
    downsample_none_dir = os.path.join("04_01_variant_analysis", "downsample_none")

    # Set output directory.
    out_root = create_dir(root_dir, "04_02_variant_analysis")
    out_downsample_dir = create_dir(out_root, f"downsample_{downsample_num}")


    bins = ["G", "L", "Y"]


    # Get all single double variants files
    single_double_variants_files = (get_files(os.path.join(downsample_dir, "04_01_03_single_double_variants", "*.parquet")) +
                                    get_files(os.path.join(downsample_none_dir, "04_01_03_single_double_variants", "*.parquet")))

    for bin in bins:
        files = [file for file in single_double_variants_files if bin in file.split("/")[-1]]
        R1_merged_df = merge_df(files, "R1")
        R2_merged_df = merge_df(files, "R2")

        R1_filtered_df = filter_df(R1_merged_df, threshold, bin, "R1")
        R2_filtered_df = filter_df(R2_merged_df, threshold, bin, "R2")


        del R1_merged_df, R2_merged_df

    # Generate input files for 05_dSortSeq_FP64.py
    abundance_filtered_files = get_files(os.path.join(out_downsample_dir, f"04_02_01_filtered_abundance_{threshold}", "*.csv"))

    for file in abundance_filtered_files:
        bin = file.split("/")[-1].split(".")[0]
        logger.info("Generating input files for 05_dSortSeq.py...")
        binned_distribution(file, bin)
        mixing_coefficient(file, bin)

    flow_csv_files = get_files(os.path.join(root_dir, "raw_flow_data", "*.csv"))
    sorting_boundaries(flow_csv_files)


    logger.info("Done!")