#!/usr/bin/env python
# coding: utf-8

# Basic diagnostics

import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
import warnings
warnings.filterwarnings('ignore')


# Import files

def import_files() -> pd.DataFrame:
    df_nuclear = pd.read_csv("data/transcripts_cellpose.csv")
    df_baysor = pd.read_csv("data/transcripts.csv")
    df_nuclear = df_nuclear[['transcript_id', 'cell_id']]
    df = pd.merge(left=df_baysor, right=df_nuclear, how="left", left_on="transcript_id", right_on='transcript_id')
    return df

df = import_files()


# Plotting functions

def diagnostic_plots(df, ax1, ax2, ax3, ax4): 
    df_cell_counts = pd.DataFrame({
        "Baysor": [df[~pd.isnull(df.cell)].cell.nunique()],
        "Nuclear\nSegmentation": [df[df.cell_id > 0].cell_id.nunique()]
    })
    sns.barplot(df_cell_counts, ax=ax1)
    ax1.set_title(f"Number of detected cells")

    # how many transcripts are assigned to a cell
    df_assigned = pd.DataFrame({
    "Baysor": (~pd.isnull(df.cell)).value_counts(),
    "Nuclear\nSegmentation": (df.cell_id > 0).value_counts()
    })
    df_assigned = df_assigned.transpose().iloc[:, [1,0]]
    df_assigned.plot.bar(stacked=True, ax=ax2, color=["#55A868", "#BD4B4F"]).legend(loc='upper right')
    ax2.set_title(f"Transcripts assigned to cell")
   
    # how many transcripts per cell
    baysor = df[~pd.isnull(df.cell)].groupby("cell").size()
    nucleus = df[df.cell_id > 0].groupby("cell_id").size()
    sns.histplot(baysor, ax=ax3, kde=True, label="Baysor", binwidth=10)
    sns.histplot(nucleus, ax=ax3, kde=True, label="Nuclear segmentation",binwidth=10)
    ax3.set_xlabel(r"Transcripts per cell")
    ax3.set_ylabel(r"cells")
    ax3.set_title(f"Transcripts per cell")
    ax3.legend()

    # how many features per cell
    baysor = df[~pd.isnull(df.cell)].groupby('cell')['gene'].nunique()
    nucleus = df[df.cell_id > 0].groupby('cell')['gene'].nunique()
    sns.histplot(baysor, ax=ax4, kde=True, label="Baysor", binwidth=5)
    sns.histplot(nucleus, ax=ax4, kde=True, label="Nuclear segmentation",binwidth=5)    
    ax4.set_xlabel(r"Features per cell")
    ax4.set_ylabel(r"cells")
    ax4.set_title("Features per cell")
    ax4.legend() 



# Make the plots

fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize = (14,5))
diagnostic_plots(df, ax1, ax2, ax3, ax4)
fig.tight_layout(h_pad=1.1)
