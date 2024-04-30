import os
from itertools import cycle  # TODO: colorlist
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statannotations.Annotator import Annotator


def concatenate_files(file_list):
    """Concatenate multiple qPCR data files into a single dataframe.

    Parameters
    ----------
    file_list : list
        List of file paths to concatenate.

    Returns
    -------
    df : pd.DataFrame
        Concatenated dataframe containing data from all input files.
    """

    df = pd.concat(
        [
            pd.read_excel(file, sheet_name="Results", skiprows=7).iloc[:-5]
            for file in file_list
            if os.path.exists(file)
        ]
    )

    if df.empty:
        warn("No data found in the input files.")

    return df


def concatenate_files_(file_dir: str):

    if not os.path.exists(file_dir):
        raise FileNotFoundError("The specified directory does not exist.")

    df = pd.concat(
        [
            pd.read_excel(file, sheet_name="Results", skiprows=7).iloc[:-5]
            for file in os.listdir(file_dir)
            if file.endswith(".xls")
        ]
    )

    if df.empty:
        warn("No data found in the input files.")

    return df


def perform_ddct_analysis(
    data_or_file: str | pd.DataFrame,
    internal_control: str,
    control_sample: str,
):
    """Perform ΔΔCt analysis on qPCR data.

    Parameters
    ----------
    data_or_file : pd.DataFrame | str
        A pandas DataFrame containing the data or path to the qPCR data file.
    internal_control : str
        Internal control gene name, like "ACTB".
    control_sample : str
        Control sample name, like "DMSO".

    Returns
    -------
    df : pd.DataFrame
        Processed dataframe.
    """

    if isinstance(data_or_file, pd.DataFrame):
        df = data_or_file
    elif isinstance(data_or_file, str) and (
        data_or_file.endswith(".xlsx") or data_or_file.endswith(".xls")
    ):
        df = pd.read_excel(data_or_file, sheet_name="Results", skiprows=7)
        df = df.iloc[:-5]
    else:
        raise ValueError(
            "Unsupported file format. Please provide a .csv or .xlsx file."
        )

    # Calculate internal control average Ct for each sample
    control_ct = (
        df[df["Target Name"] == internal_control]
        .groupby("Sample Name")["Cт"]
        .mean()
        .reset_index(drop=True)
    )
    control_ct.rename(columns={"Cт": "Internal control average Cт"}, inplace=True)

    # match internal control average Ct for each sample's genes
    df = df.merge(control_ct, on="Sample Name", how="left")

    # Calculate each sample's gene ΔCt
    df["ΔCт"] = df["Cт"] - df["Internal control average Cт"]

    # Calculate the average ΔCt for the control sample for each target gene
    control_dct = (
        df[df["Sample Name"] == control_sample]
        .groupby("Target Name")["ΔCт"]
        .mean()
        .reset_index(drop=True)
    )
    control_dct.rename(columns={"ΔCт": "Average ΔCт Control"}, inplace=True)

    df = df.merge(control_dct, on="Target Name", how="left")

    df["ΔΔCт"] = df["ΔCт"] - df["Average ΔCт Control"]

    # Calculate the relative expression level as 2^-ΔΔCt for each sample
    df["Relative Expression"] = 2 ** df["ΔΔCт"]

    return df


def plot_results_with_significance(df, internal_control, control_sample, **kwargs):
    """Plot the qPCR results with significance annotations.

    Parameters
    ----------
    df : pd.DataFrame
        Processed dataframe.
    internal_control : str
        Internal control gene name, like "ACTB".
    control_sample : str
        Control sample name, like "DMSO".
    **kwargs : dict
        sample_name : str, optional
            Name of the compound to be displayed as title.

    """

    sample_name = kwargs.get("sample_name", "Compound name")
    save_path = kwargs.get("save_path", None)

    # Exclude the internal control gene, i.e. "ACTB" will not be plotted
    plot_data = df[df["Target Name"] != internal_control].copy()

    # Change sample names into concentrations, e.g. "DMSO" -> "DMSO", "2437-50_0.1" -> "0.1 μM", to be used as hue order/legend
    concentrations: pd.Series = (
        df["Sample Name"]
        .drop_duplicates()
        .apply(lambda x: x.split("_")[-1])
        .drop_duplicates()
        .astype(float)
        .sort_values()
    )
    plot_data["Sample Name"] = plot_data["Sample Name"].map(change_sample_name)
    hue_order = concentrations.apply(lambda x: f"{x:1f} μM").to_list()
    hue_order.insert(0, "DMSO")

    # Create the bar plot using seaborn
    pre_cp = plt.get_cmap("BrBG")(np.linspace(0, 1, 100))
    if len(hue_order) == 2:
        cp = sns.color_palette([pre_cp[27], pre_cp[78]])
    elif len(hue_order) == 3:
        cp = sns.color_palette([pre_cp[27], pre_cp[68], pre_cp[78]])
    else:
        color_list = [
            pre_cp[i] for i in np.linspace(59, 99, len(hue_order)).astype(int)
        ]
        color_list.insert(0, pre_cp[27])
        cp = sns.color_palette(color_list)

    # Define gene order for plotting
    # gene_order = ['Myc', 'LDHA', 'TERT', 'ODC1', 'CDK4', 'CCND2', 'FABP5', 'G3BP1', 'SLC16A1', 'NCL', 'HK2', 'BZW2', 'FBXO32']   # GAN, FBXO32
    gene_order = ["TP53", "MDM2", "p21", "PUMA"]
    unique_genes: set = plot_data[plot_data["Sample Name"] != control_sample][
        "Target Name"
    ].to_set()
    plot_order = [gene for gene in gene_order if gene in unique_genes]

    # Initialize the matplotlib figure
    figure_width = 8 if len(plot_data) < 48 else (len(plot_data) - 1) / 5.5
    plt.figure(figsize=(figure_width, 6))
    barplot = sns.barplot(
        data=plot_data,
        x="Target Name",
        y="Relative Expression",
        hue="Sample Name",
        order=plot_order,
        hue_order=hue_order,
        palette=cp,
        errwidth=1.5,
        capsize=0.03,
        edgecolor="black",
        linewidth=1,
        saturation=1,
    )

    # set bar width and position
    for patch in barplot.patches:
        current_width = patch.get_width()
        new_width = current_width * 0.8
        diff = current_width - new_width
        patch.set_width(new_width)
        patch.set_x(patch.get_x() + diff * 0.5)

    # Add significance annotations
    # Define the box pairs to be compared
    box_pairs = [
        ((f"{target}", f"{control_sample}"), (f"{target}", f"{sample}"))
        for sample in hue_order
        if sample != control_sample
        for target in plot_order
        if target != internal_control
    ]
    annotator = Annotator(
        ax=barplot,
        pairs=box_pairs,
        data=plot_data,
        x="Target Name",
        y="Relative Expression",
        hue="Sample Name",
        order=plot_order,
        hue_order=hue_order,
    )
    annotator.configure(
        test="t-test_welch",
        text_format="star",
        line_height=0.03,
        line_width=1,
        hide_non_significant=True,
    )
    annotator.apply_and_annotate()

    plt.title(f"{sample_name}", fontsize=16, fontweight="bold", family="Arial")
    plt.xlabel("Gene name", fontsize=14, fontweight="bold", family="Arial")
    plt.ylabel("Relative mRNA level", fontsize=14, fontweight="bold", family="Arial")
    plt.xticks(fontsize=12, family="Arial")
    plt.yticks(fontsize=12, family="Arial")
    plt.legend(
        bbox_to_anchor=(1.05, 1), loc="upper left", prop={"family": "Arial", "size": 14}
    )

    # save the plot
    if save_path:
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        plt.savefig(
            os.path.join(save_path, f"{sample_name}"),
            bbox_inches="tight",
            dpi=300,
            transparent=True,
        )

    plt.tight_layout()
    plt.show()
    plt.close()


def get_sample_id(df):
    """Get sample id from the dataframe.

    Parameters
    ----------
    df : pd.DataFrame

    Returns
    -------
    sample_id : np.array
        contain unique sample id.
    """

    return pd.Series(df["Sample Name"].unique()).str.split("_").str[0].unique()


def change_sample_name(x):
    """Replace sample name with sample concentration. Concentration should be marked after sample name with '_'.

    Parameters
    ----------
    x : pd.Series.value

    Returns
    -------
    y : str
    """

    return "DMSO" if x.split("_")[0] == "DMSO" else x.split("_")[1] + " μM"


def change_float_to_str(x):
    """Change float to string with 2 decimal places.

    Parameters
    ----------
    x : float

    Returns
    -------
    y : str
    """

    # if float is a integer, return as integer
    return str(int(x)) if x == int(x) else str(round(x, 2))
