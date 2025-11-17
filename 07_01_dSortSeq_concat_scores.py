import allel
import glob
import pandas as pd
from collections import Counter
from datetime import datetime
import argparse
import os
import logging
import cudf
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import dask.dataframe as dd
import plotly.express as px
from scipy.special import softmax
import plotly.graph_objects as go
from tensorboard.data.proto.data_provider_pb2 import Downsample


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


def get_downsample_files(start, end):
    """
    Get dSortSeq score csv files from downsample directories.
    :param start: starting downsample folder number
    :param end: ending downsample folder number
    :return: a list of dSortSeq score csv files
    """
    downsample_files = []
    logger.info(f"Getting downsample files from {start} to {end}...")
    for i in range(start, end+1):
        downsample_dir = os.path.join(root_dir, "06_01_dSortSeq_results", f"downsample_{i}")
        downsample_files += get_files(f"{downsample_dir}/*.csv")

    return downsample_files


def concat_scores(files):
    """
    Concatenate dSortSeq score csv files.
    :param files: list of files
    :return: merged dataframe with mean and sd columns
    """

    dark_files = [file for file in files if "dark" in file]
    light_files = [file for file in files if "light" in file]
    ymScarlet_files = [file for file in files if "ymScarlet" in file]

    dark_df = csv_to_df(dark_files)
    light_df = csv_to_df(light_files)
    ymScarlet_df = csv_to_df(ymScarlet_files)

    dfs = {"dark_df": dark_df, "light_df": light_df, "ymScarlet_df": ymScarlet_df}

    activation_df = cal_activation_df(dfs)

    dfs["activation_df"] = activation_df


    return dfs


def cal_activation_df(dfs):

    temp_df = pd.merge(dfs["dark_df"], dfs["light_df"], on="variant", how="inner", suffixes=("_dark", "_light"))

    activation_df = temp_df.copy()
    activation_df["mean"] = activation_df["mean_light"] / activation_df["mean_dark"]
    activation_df["variant_type"] = activation_df["variant_type_dark"]

    activation_df = activation_df[["variant", "mean", "variant_type"]]
    return activation_df


def csv_to_df(files):
    df = pd.DataFrame([], columns=['variant'])
    for i in range(len(files)):
        temp_df = pd.read_csv(files[i], usecols=['variant', 'mean'])
        temp_df.columns = ["variant", f"mean{i+1}"]
        # Outer join dataframes.
        df = pd.merge(df, temp_df, on='variant', how='outer')
    # Create a new column with the average score.
    df['mean'] = df.iloc[:, 1:].mean(axis=1)

    # Remove outliers using IQR method.

    # Calculate Q1, Q3, and IQR for each variant.
    mean_cols = [f"mean{i+1}" for i in range(len(files))]
    stats_df = pd.DataFrame([], columns=['variant'])
    stats_df['variant'] = df['variant']
    stats_df['Q1'] = df[mean_cols].quantile(0.25, axis=1)
    stats_df['Q3'] = df[mean_cols].quantile(0.75, axis=1)
    stats_df['IQR'] = stats_df['Q3'] - stats_df['Q1']
    stats_df['lower_bound'] = stats_df['Q1'] - 1.5 * stats_df['IQR']
    stats_df['upper_bound'] = stats_df['Q3'] + 1.5 * stats_df['IQR']

    # Remove outliers in df mean_cols using IQR method (Q1 - 1.5 * IQR, Q3 + 1.5 * IQR)
    df_filtered = pd.merge(df, stats_df, on='variant', how='left')

    for col in mean_cols:
        mask = (df_filtered[col] < df_filtered['lower_bound']) | (df_filtered[col] > df_filtered['upper_bound'])

        # Replace outliers with NaN
        df_filtered[col] = df_filtered[col].mask(mask, np.nan)

    # Recalculate mean and sd
    df_filtered['mean'] = df_filtered[mean_cols].mean(axis=1)

    # Get variant type

    def get_variant_type(variant):
        if "*" in variant:
            return "Nonsense"
        elif variant[0] == variant[-1]:
            return "Synonymous"
        else:
            return "Missense"

    df_filtered["variant_type"] = df_filtered["variant"].apply(get_variant_type)

    return df_filtered[['variant', 'mean', 'variant_type']]


def export_dfs(dfs):

    for condition in ["dark", "light", "ymScarlet", "activation"]:
        df = dfs[f"{condition}_df"]
        df.to_csv(os.path.join(output_dir, f"{condition}_df.csv"), index=False)

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
            title=f"{condition} score",
            xaxis_title="Position",
            yaxis_title="Amino acid"
        )


        fig.write_html(os.path.join(output_dir, f"{condition}_heatmap.html"))
        # fig.show()


    return None


def hist_plots(dfs):

    colour_max = {}
    colour_min = {}

    for condition in ["dark", "light", "ymScarlet", "activation"]:
        colour_max[condition] = dfs[f"{condition}_df"]["mean"].quantile(0.998)
        colour_min[condition] = dfs[f"{condition}_df"]["mean"].quantile(0.005)

    colour_max["dark"] = colour_max["light"]
    colour_min["dark"] = colour_min["light"]

    for condition in ["dark", "light", "ymScarlet", "activation"]:

        df = dfs[f"{condition}_df"]

        fig = px.histogram(
            df,
            x="mean",
            color="variant_type",
            barmode="overlay",
            hover_data=["variant"]
        )

        fig.update_layout(
            title=f"{condition} score",
            xaxis_title="Mean",
            yaxis_title="Count",
            xaxis_range=[colour_min[condition], colour_max[condition]]
        )

        fig.show()

    return None

def scatter_plots(dfs):
    activation_df = dfs["activation_df"][["variant", "mean"]].rename(columns={"mean": "activation"})
    ymScarlet_df = dfs["ymScarlet_df"][["variant", "mean", "variant_type"]].rename(columns={"mean": "ymScarlet"})
    data = pd.merge(activation_df, ymScarlet_df, on="variant", how="inner")

    x_max = data["activation"].quantile(0.998)
    x_min = data["activation"].quantile(0.005)
    y_max = data["ymScarlet"].quantile(0.998)
    y_min = data["ymScarlet"].quantile(0.005)

    fig = px.scatter()

    fig.add_scatter(
        x=data[data["variant_type"] == "Missense"]["activation"],
        y=data[data["variant_type"] == "Missense"]["ymScarlet"],
        mode="markers",
        name="Missense",
        marker=dict(color="black"),
        hovertext=data[data["variant_type"] == "Missense"]["variant"]
    )

    fig.add_scatter(
        x=data[data["variant_type"] == "Nonsense"]["activation"],
        y=data[data["variant_type"] == "Nonsense"]["ymScarlet"],
        mode="markers",
        name="Nonsense",
        marker=dict(color="red"),
        hovertext=data[data["variant_type"] == "Nonsense"]["variant"]
    )

    fig.add_scatter(
        x=data[data["variant_type"] == "Synonymous"]["activation"],
        y=data[data["variant_type"] == "Synonymous"]["ymScarlet"],
        mode="markers",
        name="Synonymous",
        marker=dict(color="blue"),
        hovertext=data[data["variant_type"] == "Synonymous"]["variant"]
    )

    fig.update_layout(
        template="simple_white",
        xaxis_title="Activation",
        yaxis_title="Per cell protein abundance",
        xaxis_range=[x_min, x_max],
        yaxis_range=[y_min, y_max],
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.01,
            xanchor="right",
            x=1,
            font=dict(size=16)
        ),
        autosize=False,
        width=700,
        height=600
    )

    fig.update_xaxes(linewidth=2.5,
                     linecolor="black",
                     tickwidth=2.5,
                     title_font=dict(size=18),
                     tickfont=dict(size=16)
                     )
    fig.update_yaxes(linewidth=2.5,
                     linecolor="black",
                     tickwidth=2.5,
                     title_font=dict(size=18),
                     tickfont=dict(size=16)
                     )

    fig.update_traces(
        hovertemplate="Variant: %{hovertext}<br>Activation: %{x:.2f}<br>Per cell protein abundance: %{y:.2f}<extra></extra>")

    fig.show()



if __name__ == '__main__':

    # Create log file.
    logger = create_log(False)

    # Get root directory.
    root_dir = os.getcwd()
    logger.info(f"Root directory: {root_dir}")

    # Get arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--start", help="Downsample folder number (starting)", type=int)
    parser.add_argument("-M", "--end", help="Downsample folder number (ending)", type=int)
    args = parser.parse_args()

    downsample_num_start = args.start
    downsample_num_end = args.end


    # Set output root directory
    output_root_dir = create_dir(root_dir, "07_01_concat_dSortSeq_scores")
    output_dir = create_dir(output_root_dir, f"downsample_{downsample_num_start}_{downsample_num_end}")

    # Get dSortSeq score csv files from downsample directories.
    downsample_files = get_downsample_files(downsample_num_start, downsample_num_end)

    # Concatenate dSortSeq score csv files
    dfs = concat_scores(downsample_files)


    # Export dataframes
    export_dfs(dfs)

    # Create heat maps
    heat_maps(dfs)


    logger.info("Done!")