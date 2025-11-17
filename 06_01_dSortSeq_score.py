import glob
import pandas as pd
from datetime import datetime
import argparse
import os
import logging
import numpy as np
from scipy.special import softmax
import plotly.graph_objects as go


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


def get_results():

    conditions = ["dark", "light", "ymScarlet"]
    dfs = {}

    for condition in conditions:
        variants = pd.read_csv(os.path.join(downsample_inputs_dir, f"{condition}_binned_distribution.csv")).iloc[:,[0]]
        variants.columns = ["variant"]
        # Extract single variants from the variant column. Single variants = len(variant.split("_")) == 1
        single_variants = variants[variants["variant"].apply(lambda x: len(x.split("_")) == 1)]

        results_dir = os.path.join(downsample_outputs_dir, condition)

        mu = pd.read_pickle(os.path.join(results_dir, "mu.pickle"))
        sigma = pd.read_pickle(os.path.join(results_dir, "sigma.pickle"))
        lamb = pd.read_pickle(os.path.join(results_dir, "lamb.pickle"))
        lamb = softmax(lamb, axis=1)

        m = mu * np.log(10 ** lamb)
        V = (sigma * np.log(10 ** lamb)) ** 2
        mean = np.exp(m + V / 2)

        Mean = mean.prod(axis=1)
        Mean = pd.DataFrame(Mean, columns=["mean"])


        # Merge dataframes based on index
        dfs[f"{condition}_df"] = pd.concat([variants, Mean], axis=1)
        # Only keep single variants
        dfs[f"{condition}_df"] = dfs[f"{condition}_df"].merge(single_variants, on="variant", how="inner")


    # Get the activation dataframe
    dfs["activation_df"] = cal_activation_df(dfs)

    def get_variant_type(variant):
        if "*" in variant:
            return "Nonsense"
        elif variant[0] == variant[-1]:
            return "Synonymous"
        else:
            return "Missense"

    for df in dfs.values():
        df["variant_type"] = df["variant"].apply(get_variant_type)

    return dfs


def cal_activation_df(dfs):

    temp_df = pd.merge(dfs["dark_df"], dfs["light_df"], on="variant", how="inner", suffixes=("_dark", "_light"))

    activation_df = temp_df.copy()
    activation_df["mean"] = activation_df["mean_light"] / activation_df["mean_dark"]

    activation_df = activation_df[["variant", "mean"]]
    return activation_df


def heat_maps(dfs):

    dfs = prep_df_for_heat_maps(dfs)

    colour_max = {}
    colour_min = {}

    for condition in ["dark", "light", "ymScarlet", "activation"]:
        colour_max[condition] = dfs[f"{condition}_df"]["mean"].quantile(0.998)
        colour_min[condition] = dfs[f"{condition}_df"]["mean"].quantile(0.005)

    colour_max["dark"] = colour_max["light"]
    colour_min["dark"] = colour_min["light"]

    # Create a heat map for each df in dfs. The row is the position of the variant and the column is the MUT amino acid.
    for condition in ["dark", "light", "ymScarlet", "activation"]:

        df = dfs[f"{condition}_df"]


        # Create and save heat map for mean. X axis is the position of the variant and Y axis is the amino acid. The
        # colour is the mean value. Show df['variant'] as hover data.
        fig = go.Figure(data=go.Heatmap(x=df["position"],
                                        y=df["MUT"],
                                        z=df["mean"],
                                        zmin=colour_min[condition],
                                        zmax=colour_max[condition],
                                        hovertext=df["variant"],
                                        hovertemplate="Variant: %{hovertext}<br>Mean: %{z:.2f}<extra></extra>"
                                        )
                         )

        # Order y axis by volume of amino acids
        fig.update_yaxes(categoryorder="array",
                         categoryarray=["*",
                                        "G", "A", "S", "C", "D",
                                        "P", "N", "T", "E", "V",
                                        "Q", "H", "M", "I", "L",
                                        "K", "R", "F", "Y", "W"])


        fig.update_layout(
            title=f"Downsample {downsample_num} {condition} score",
            xaxis_title="Position",
            yaxis_title="Amino acid"
        )



        fig.write_html(os.path.join(output_dir, f"{condition}_heatmap.html"))
        # scatter_fig.show()


    return None


def prep_df_for_heat_maps(dfs):

    for condition in ["dark", "light", "ymScarlet", "activation"]:
        df = dfs[f"{condition}_df"]
        # Get position of variant by getting the digit in the variant string
        df["position"] = df["variant"].str.extract(r'(\d+)').astype(int)
        # Get WT which is the first character in the variant string
        df["WT"] = df["variant"].str[0]
        # Get MUT which is the last character in the variant string
        df["MUT"] = df["variant"].str[-1]
        # Sort by position and reset index
        df = df.sort_values("position").reset_index(drop=True)

        # Change variant column object to dictionary with name, position, WT, and MUT
        df = df[["variant", "position", "WT", "MUT", "mean", "variant_type"]]

    return dfs



def export_dfs(dfs):

    for condition in ["dark", "light", "ymScarlet", "activation"]:
        df = dfs[f"{condition}_df"]
        df.to_csv(os.path.join(output_dir, f"{condition}_df.csv"), index=False)

    return None


if __name__ == '__main__':

    # Create log file.
    logger = create_log(False)

    # Get root directory.
    root_dir = os.getcwd()
    logger.info(f"Root directory: {root_dir}")

    # Get arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--number", help="Downsample folder number", type=int)
    args = parser.parse_args()

    downsample_num = args.number

    downsample_inputs_dir = os.path.join(root_dir, "04_01_dSortSeq_inputs", f"downsample_{downsample_num}")
    downsample_outputs_dir = os.path.join(root_dir, "05_01_dSortSeq_outputs", f"downsample_{downsample_num}")

    output_root_dir = create_dir(root_dir, "06_01_dSortSeq_results")
    output_dir = create_dir(output_root_dir, f"downsample_{downsample_num}")

    dfs = get_results()

    export_dfs(dfs)

    heat_maps(dfs)

    logger.info("Done!")