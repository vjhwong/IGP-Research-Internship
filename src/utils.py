import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import kaplanmeier as km


"""Functions used in sia_survivalplot.ipynb"""


def plot_expression_histograms(expression_df) -> None:
    _, axs = plt.subplots(2, 2, figsize=(10, 10))
    for i, gene in enumerate(expression_df.columns):
        sns.histplot(expression_df[gene], bins=50, ax=axs[i // 2, i % 2])
        axs[i // 2, i % 2].set_xlabel("Values")
        axs[i // 2, i % 2].set_ylabel("Frequency")
        axs[i // 2, i % 2].set_title(f"Histogram of {gene} expression values")
    plt.show()


def format_data(expression: pd.DataFrame, C1QX: str) -> pd.DataFrame:
    """Formats the data for the SIA calculation and categorization

    Args:
        expression (pd.Dataframe): Contains the expression data
        genes (list[str]): List of genes to be used for the SIA calculation

    Returns:
        pd.Dataframe: Dataframe with the SIA and SIA-category {C1QX} columns
    """
    expression[f"SIA-category {C1QX}"] = pd.qcut(
        expression[f"SIA {C1QX}"], 3, labels=["Low SIA", "Intermediate SIA", "High SIA"]
    )
    return expression.dropna()


def get_log_rank_p(expression: pd.DataFrame, C1QX: str) -> tuple[float, float]:
    """Calculates the log-rank p-value for the SIA categories

    Args:
        expression (pd.DataFrame): Dataframe containing the SIA and SIA-category columns

    Returns:
        tuple[float, float]: Log-rank p-values for high vs intermediate and intermediate vs low
    """
    # Partition data
    expression_high_intermediate = expression[
        (
            (expression[f"SIA-category {C1QX}"] == "High SIA")
            | (expression[f"SIA-category {C1QX}"] == "Intermediate SIA")
        )
    ]
    expression_intermediate_low = expression[
        (
            (expression[f"SIA-category {C1QX}"] == "Intermediate SIA")
            | (expression[f"SIA-category {C1QX}"] == "Low SIA")
        )
    ]

    # Calculate log-rank p-value for high vs intermediate
    time_high_intermediate = expression_high_intermediate["os.delay"]
    event_high_intermediate = expression_high_intermediate["os.event"]
    group_high_intermediate = expression_high_intermediate[f"SIA-category {C1QX}"]
    results_high_intermediate = km.fit(
        time_high_intermediate, event_high_intermediate, group_high_intermediate
    )
    log_rank_p_high_intermediate = results_high_intermediate["logrank_P"]

    # Calculate log-rank p-value for intermediate vs low
    time_intermediate_low = expression_intermediate_low["os.delay"]
    event_intermediate_low = expression_intermediate_low["os.event"]
    group_intermediate_low = expression_intermediate_low[f"SIA-category {C1QX}"]
    results_intermediate_low = km.fit(
        time_intermediate_low, event_intermediate_low, group_intermediate_low
    )
    log_rank_p_intermediate_low = results_intermediate_low["logrank_P"]

    return log_rank_p_high_intermediate, log_rank_p_intermediate_low


def plot_survival_curve(
    expression: pd.DataFrame,
    log_rank_p_high_intermediate: float,
    log_rank_p_intermediate_low: float,
    C1QX: str,
) -> None:
    """Plots the survival curve for the SIA categories

    Args:
        expression (pd.DataFrame): Dataframe containing the SIA and SIA-category columns
    """
    high_intermediate_label = (
        f"High vs Intermediate (p = {log_rank_p_high_intermediate:.3f})"
    )
    intermediate_low_label = (
        f"Intermediate vs Low (p = {log_rank_p_intermediate_low:.3f})"
    )
    plot_title = (
        f"Survival curve for SIA categories "
        + f"({high_intermediate_label}, {intermediate_low_label})"
    )
    time = expression["os.delay"]
    event = expression["os.event"]
    group = expression[f"SIA-category {C1QX}"]
    results = km.fit(time, event, group)
    km.plot(
        results,
        title=plot_title,
    )


"""Function used in sia_stripplot.ipynb"""


def process_plot_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Process the data for plotting. This includes combining Pt27A and Pt27B into one row, removing Pt20 and renaming the values in response_ordinal column.

    Args:
        df (pd.DataFrame): Dataframe to process for plotting.

    Returns:
        pd.DataFrame: Processed dataframe for plotting.
    """
    # Rename the values in response_ordinal column
    response_mapper = {
        "Progressive Disease": "non-responder",
        "Partial Response": "partial response",
        "Complete Response": "responder",
    }
    df["response"] = df["response_ordinal"].apply(lambda x: response_mapper[x])

    # Columns used for plotting
    plot_columns = [
        "SIAC1QA",
        "SIAC1QB",
        "SIAC1QC",
        "SIA celltype",
        "response",
        "ID_REF",
    ]
    plot_data = df[plot_columns]

    # Combine Pt27A and Pt27B into one row
    row_index_Pt27A = df[df["ID_REF"] == "Pt27A"].index[0]
    row_index_Pt27B = df[df["ID_REF"] == "Pt27B"].index[0]
    new_row = pd.DataFrame(
        {
            # Take the mean of the values for Pt27A and Pt27B
            "SIAC1QA": np.array(
                [
                    plot_data.loc[row_index_Pt27A, "SIAC1QA"],
                    plot_data.loc[row_index_Pt27B, "SIAC1QA"],
                ]
            ).mean(),
            "SIAC1QB": np.array(
                [
                    plot_data.loc[row_index_Pt27A, "SIAC1QB"],
                    plot_data.loc[row_index_Pt27B, "SIAC1QB"],
                ]
            ).mean(),
            "SIAC1QC": np.array(
                [
                    plot_data.loc[row_index_Pt27A, "SIAC1QC"],
                    plot_data.loc[row_index_Pt27B, "SIAC1QC"],
                ]
            ).mean(),
            "SIA celltype": np.array(
                [
                    plot_data.loc[row_index_Pt27A, "SIA celltype"],
                    plot_data.loc[row_index_Pt27B, "SIA celltype"],
                ]
            ).mean(),
            "response": plot_data.loc[row_index_Pt27A, "response"],
            "ID_REF": "Pt27A_Pt27B",
        },
        index=[0],
    )
    # Remove Pt27A, Pt27B and outlier Pt20
    plot_data = plot_data.drop(
        plot_data[plot_data["ID_REF"].isin(["Pt27A", "Pt27B", "Pt20"])].index
    )

    # Add the new row
    plot_data = pd.concat([plot_data, new_row], ignore_index=True)
    return plot_data


def find_min_responder(df: pd.DataFrame, column: str) -> float:
    """
    Find the minimum value for the responders in the given column.

    Args:
        df (pd.Series): Dataframe containing the data.
        column (str): Column to find the minimum value for.

    Returns:
        float: Minimum value for the responders in the given column.
    """
    return df[df["response"] == "responder"][column].min()


def count_less_than_min_responder(
    df: pd.DataFrame, column: str, minimum_responder: float
) -> dict[str:int]:
    """
    Count the number of non-responders and partial responders with values less than the minimum responder.

    Args:
        df (pd.DataFrame): Dataframe containing the data.
        column (str): Column to count the number of non-responders and partial responders for.
        minimum_responder (float): Minimum value for the responders in the given column.

    Returns:
        dict: Dictionary containing the counts for non-responders and partial responders with values less than the minimum responder.
    """
    count = {"non_respoder": 0, "partial_responder": 0}

    # Select the rows where the response is non_responder and partial_responder
    count["non_responder"] = df[
        (df["response"] == "non-responder") & (df[column] < minimum_responder)
    ].shape[0]
    count["partial_responder"] = df[
        (df["response"] == "partial response") & (df[column] < minimum_responder)
    ].shape[0]
    return count


def generate_stripplot(plot_data: pd.DataFrame) -> None:
    """
    Generate a stripplot for each column in the plot_data DataFrame.

    Args:
        plot_data (pd.DataFrame): The DataFrame to plot.
    """
    fig, axes = plt.subplots(1, len(plot_data.columns) - 2, figsize=(20, 5))
    n_non_responders = plot_data[plot_data["response"] == "non-responder"].shape[0]
    n_partial_responders = plot_data[plot_data["response"] == "partial response"].shape[
        0
    ]

    for i, column in enumerate(plot_data.columns[:-2]):
        min_responder = find_min_responder(plot_data, column)
        count = count_less_than_min_responder(plot_data, column, min_responder)
        sns.stripplot(
            x="response",
            y=column,
            data=plot_data,
            jitter=True,
            palette="Set1",
            hue="response",
            ax=axes[i],
            size=10,
            alpha=0.7,
            edgecolor="black",
            linewidth=1,
        )
        axes[i].axhline(y=min_responder, color="g", linestyle="--")
        axes[i].set_xlabel("")
        axes[i].set_ylabel("SIA")
        if column.startswith("SIAC1Q"):
            axes[i].set_ylim(0, 0.1)
        test_str = "\n".join(
            (
                f"NR less than MR: {count['non_responder']}/{n_non_responders}",
                f"PR less than MR: {count['partial_responder']}/{n_partial_responders}",
            )
        )
        props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
        axes[i].text(
            0.05,
            0.95,
            test_str,
            transform=axes[i].transAxes,
            fontsize=10,
            verticalalignment="top",
            bbox=props,
        )
        axes[i].text(
            1.25,
            min_responder * 0.6,
            f"Min. Responder = {min_responder:.3f}",
            fontsize=9,
            color="g",
        )
        axes[i].set_title(column)
    fig.suptitle("Melanoma (bulk RNA)")
    # Adjust layout
    plt.tight_layout()

    # Show plot
    plt.show()


def plot_gene_celltype_correlation(
    data_qs: pd.DataFrame,
    deconvolved_data: pd.DataFrame,
    genes: list[str],
    cell_type: str,
) -> None:
    """
    Create a heatmap of the correlation between gene expression and cell type.

    Args:
        data_qs (pd.DataFrame): Gene expression data
        deconvolved_data (pd.DataFrame): Deconvolved cell type data
        genes (list[str]): List of genes
        cell_type (str): Cell type to correlate with
    """
    # Select the genes from the gene expression data
    df = data_qs[data_qs["Symbol"].isin(genes)]
    df.set_index("Symbol", inplace=True)
    df = df.T
    df.reset_index(inplace=True, names=None)
    df.rename(columns={"index": "Mixture"}, inplace=True)

    # Merge the gene expression data with the deconvolved cell type data
    df = df.merge(deconvolved_data[["Mixture", cell_type]], on="Mixture")
    df.drop(columns=["Mixture"], inplace=True)

    # Create the heatmap
    sns.heatmap(df.corr(), annot=True, cmap="coolwarm")
    plt.title("Gene Expression and Celltype Heatmap")
    plt.show()
