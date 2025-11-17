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
        downsample_dir = os.path.join(root_dir, "06_02_dSortSeq_results", f"downsample_{i}")
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
    activation_df["R1"] = activation_df["R1_light"] / activation_df["R1_dark"]
    activation_df["R2"] = activation_df["R2_light"] / activation_df["R2_dark"]
    activation_df["variant_type"] = activation_df["variant_type_dark"]

    activation_df = activation_df[["variant", "R1", "R2", "variant_type"]]
    return activation_df


def csv_to_df(files):
    replicate_dfs = {}
    for replicate in ["R1", "R2"]:
        df = pd.DataFrame([], columns=['variant'])
        for i in range(len(files)):
            temp_df = pd.read_csv(files[i], usecols=['variant', replicate])
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
        df_filtered = pd.merge(df, stats_df, on='variant', how='inner')

        for col in mean_cols:
            mask = (df_filtered[col] < df_filtered['lower_bound']) | (df_filtered[col] > df_filtered['upper_bound'])

            # Replace outliers with NaN
            df_filtered[col] = df_filtered[col].mask(mask, np.nan)

        # Recalculate mean and sd
        df_filtered[replicate] = df_filtered[mean_cols].mean(axis=1)

        replicate_dfs[replicate] = df_filtered[['variant', replicate]].copy()

    df = replicate_dfs["R1"].merge(replicate_dfs["R2"], on='variant', how='inner').dropna()

    # Get variant type
    def get_variant_type(variant):
        if "*" in variant:
            return "Nonsense"
        elif variant[0] == variant[-1]:
            return "Synonymous"
        else:
            return "Missense"

    df["variant_type"] = df["variant"].apply(get_variant_type)

    return df[['variant', 'R1', 'R2', 'variant_type']]


def export_dfs(dfs):

    for condition in ["dark", "light", "ymScarlet", "activation"]:
        df = dfs[f"{condition}_df"]
        df.to_csv(os.path.join(output_dir, f"{condition}_df.csv"), index=False)

    return None


def scatter_plots(dfs):

    for condition in ["dark", "light", "ymScarlet", "activation"]:


        fig = go.Figure(data=go.Scatter(
            x=dfs[f"{condition}_df"]["R1"],
            y=dfs[f"{condition}_df"]["R2"],
            mode='markers',
            marker=dict(
                size=10,
                color=dfs[f"{condition}_df"]["variant_type"].map({
                    "Nonsense": "red",
                    "Synonymous": "blue",
                    "Missense": "green"
                }),
                opacity=0.7,
                line=dict(width=1, color='black')
            ),
            text=dfs[f"{condition}_df"]["variant"],
            hovertemplate="<b>Variant:</b> %{text}<br>" +
                          "<b>Replicate 1:</b> %{x:.2f}<br>" +
                          "<b>Replicate 2:</b> %{y:.2f}<br>" +
                          "<b>Variant Type:</b> %{marker.color}<extra></extra>"
        ))



        fig.update_layout(
            template="simple_white",
            title=f"{condition.capitalize()}",
            xaxis_title="Replicate 1",
            yaxis_title="Replicate 2",
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
                         tickfont=dict(size=16),
                         type="log",
                         range=[np.log10(dfs[f"{condition}_df"]["R1"].quantile(0.05)), np.log10(dfs[f"{condition}_df"]["R1"].quantile(0.997))]
                         )
        fig.update_yaxes(linewidth=2.5,
                         linecolor="black",
                         tickwidth=2.5,
                         title_font=dict(size=18),
                         tickfont=dict(size=16),
                         type="log",
                         range=[np.log10(dfs[f"{condition}_df"]["R1"].quantile(0.05)), np.log10(dfs[f"{condition}_df"]["R2"].quantile(0.997))]
                         )


        fig.show()
    return None


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
    output_root_dir = create_dir(root_dir, "07_02_concat_dSortSeq_scores")
    output_dir = create_dir(output_root_dir, f"downsample_{downsample_num_start}_{downsample_num_end}")

    # Get dSortSeq score csv files from downsample directories.
    downsample_files = get_downsample_files(downsample_num_start, downsample_num_end)

    # Concatenate dSortSeq score csv files
    dfs = concat_scores(downsample_files)


    # Export dataframes
    export_dfs(dfs)

    # Create scatter plots
    scatter_plots(dfs)


    logger.info("Done!")