#!/usr/bin/env python
# coding: utf-8

from . import logger, log_step

# # Load libraries
    
import cooler
from cooler import binnify, read_chromsizes
import cooltools
import bioframe
import bioframe.sandbox.gtf_io

import pandas as pd
import numpy as np
import scipy
import statsmodels.api as sm

import tqdm
from functools import partial
import sys
tqdm.tqdm = partial(tqdm.tqdm, file=sys.stdout, colour="GREEN")

import os
import glob

# Visualization settings
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Advanced plotting settings:
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

from datashader.mpl_ext import dsshow
import datashader as ds

import ultraplot

# Adjust images properties:
ultraplot.rc["figure.facecolor"] = "white"
ultraplot.rc.update(
    linewidth=1,
    fontsize=14,
    color="black",
    suptitlecolor="black",
    titleloc="upper center",
    titlecolor="dark blue",
    titleborder=False,
)
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["grid.alpha"] = 0

# Parameter for saving text as text and not vectorized objects:
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


import fontanka
import cooltools # Loading Hi-C cooler data
import bbi # Loading epigenetics bigwig data

# Handling variables in a way that do not allow the sections to interfere:
import inspect
import types
import gc

_section_baseline = None  # will store names present at section start
def section_start():
    """
    Mark the beginning of a code section.
    """
    global _section_baseline

    # Get caller's global namespace (typically the module / notebook globals)
    frame = inspect.currentframe().f_back
    g = frame.f_globals

    _section_baseline = set(g.keys())
    print("✅ Section started")


def section_flush():
    """
    Remove all variables created since the last section_start(),
    except:
      - system/dunder variables (__name__, etc.)
      - imported modules
    """
    global _section_baseline
    # if _section_baseline is None:
    #     return  # nothing to do

    frame = inspect.currentframe().f_back
    g = frame.f_globals

    baseline = _section_baseline

    for name in list(g.keys()):
        # keep variables that already existed at section_start()
        if name in baseline:
            continue

        # keep dunder / system / internal variables
        if name.startswith("__") or name.startswith("_"):
            continue

        obj = g[name]

        # keep imported modules (libraries)
        if isinstance(obj, types.ModuleType):
            continue

        # otherwise: delete newly created variable
        del g[name]

    gc.collect()
    print("✅ Section cleaned")


# Checking the versions
# List and check all the dependencies:
import importlib.metadata as md
from pathlib import Path
from typing import Optional
import platform
import re
try:
    from packaging.version import Version
except ImportError:
    Version = None  # falls back to string compare

_SPEC_RE = re.compile(r"""
    ^\s*
    (?P<name>[A-Za-z0-9_.\-]+)
    \s*
    (?P<op>==|>=|<=|>|<)?
    \s*
    (?P<ver>[^\s#]+)?
""", re.VERBOSE)

def check_dep_versions(dep_file: Optional[Path] = None):
    """
    Read specs from dep_file (one per line, e.g. 'numpy>=1.26.4' or 'pandas <=2.2.2')
    and return {name: {installed, expected, op, status}}.
    """
    dep_path = dep_file or Path(__file__).with_name("requirements.txt")
    if not dep_path.exists():
        return {"error": f"dependency file not found: {dep_path}"}

    out = {}

    # Add Python version
    python_version = platform.python_version()
    python_version_expected = "python==3.12"
    m = _SPEC_RE.match(python_version_expected)
    name = m.group("name")
    op = m.group("op") or None
    expected = m.group("ver") if m.group("ver") else None
    out['Python'] = {"installed": python_version, "expected": expected, "op": op, "status": "installed"}

    # Add versions of all packages:
    for raw in dep_path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        m = _SPEC_RE.match(line)
        if not m:
            continue
        name = m.group("name")
        op = m.group("op") or None
        expected = m.group("ver") if m.group("ver") else None

        try:
            installed = md.version(name)
        except md.PackageNotFoundError:
            out[name] = {"installed": "not installed", "expected": expected, "op": op, "status": "missing"}
            continue

        status = "installed"
        if expected and op and Version:
            iv, ev = Version(installed), Version(expected)
            if (
                (op == "==" and iv == ev)
                or (op == ">=" and iv >= ev)
                or (op == "<=" and iv <= ev)
                or (op == ">" and iv > ev)
                or (op == "<" and iv < ev)
            ):
                status = "ok"
            else:
                status = "mismatch"
        elif expected and op:  # fallback without packaging
            status = "ok" if installed == expected else "mismatch"
        out[name] = {"installed": installed, "expected": expected, "op": op, "status": status}
    return out

# TODO: transition to more transparent code, such as in cooltools examples: 
# from packaging import version
# if version.parse(cooltools.__version__) < version.parse('0.5.2'):
#     raise AssertionError("tutorials rely on cooltools version 0.5.2 or higher,"+
#                          "please check your cooltools version and update to the latest")

    
# # Set up saving mode (no saving by default), but we create a temporal folder ./test
log_step("Default saving mode set to do_save False")
do_save = False

log_step('Making sure output folder is defined and exists in the location where the library is imported: ./test')
output_folder = "./test"
if not os.path.isdir(output_folder):
    os.mkdir(output_folder)


# Define color scheme:

# Profiles wild-type
color_fountains_dome = "#C7A33B"  # yellow
color_fountains_late = "#A66D2D"  # sienna
color_insulation_early = "#A66D2D"  # rose
color_insulation_late = "#782C2C"  # dark-red

color_IZ = "#2a2b77" # blue

# Profiles knockouts:
color_fount_down = "#892EC0"
color_fount_same = "gray"
color_fount_up = "#369A2F"
color_WT = "gray"
color_MZTriple = "#F9A350"

# Colors of the developmental stages:
color_sperm = "#84B6E2"
color_2_75 = "#B0B1B1"
color_5_3 = "#FDBF63"
color_11 = "#F58454"
color_25 = "#D15951"


# # Load universal functions and libraries

# Function to annotate the boxplot by Mann-Whitney U-test
from itertools import combinations
def annotate_boxplot(
    data,
    x,
    y,
    order,
    ax,
    alt="min",
    shift_asterisks=False,
    max_level=None,
    max_asterisks=4,
    print_ns=False,
    plot_hlines=False,
):
    """Run pairwise Mann–Whitney tests across groups and draw significance marks on a boxplot."""

    groups = order  # df_selected['type'].unique()
    p_values = {}

    print('Logging: group1, group2, p-value, alternative type, p_value less, p_value greater')
    for group1, group2 in combinations(groups, 2):

        p_value_full = None
        if alt == "min":
            p_value1 = scipy.stats.mannwhitneyu(
                data.query(f'{x}=="{group1}"').loc[:, y].dropna().values.astype(float),
                data.query(f'{x}=="{group2}"').loc[:, y].dropna().values.astype(float),
                alternative="less",
            ).pvalue
            p_value2 = scipy.stats.mannwhitneyu(
                data.query(f'{x}=="{group1}"').loc[:, y].dropna().values.astype(float),
                data.query(f'{x}=="{group2}"').loc[:, y].dropna().values.astype(float),
                alternative="greater",
            ).pvalue
            p_value = min([p_value1, p_value2])
            p_value_full = [p_value1, p_value2]
        else:
            p_value = scipy.stats.mannwhitneyu(
                data.query(f'{x}=="{group1}"').loc[:, y].dropna().values.astype(float),
                data.query(f'{x}=="{group2}"').loc[:, y].dropna().values.astype(float),
                alternative=alt,
            ).pvalue
        p_values[(group1, group2)] = p_value
        print(group1, group2, p_value, alt, *p_value_full)

    y_val = ax.get_ylim()[1]  # np.nanpercentile(df_selected[y], 99.9)
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    print('Summary:', x, y, p_values, alt, y_val)

    alpha = 0.05
    significance_levels = [0.0001, 0.001, 0.01, 0.05]
    symbols = ["****", "***", "**", "*"]

    assert max_asterisks <= 4, "Max 4 asterisks supported"
    if max_asterisks < 4:
        symbols = [x[:max_asterisks] for x in symbols]

    for i, ((group1, group2), p_value) in enumerate(p_values.items()):

        x1, x2 = groups.index(group1), groups.index(group2)
        if max_level is not None:
            if abs(x1 - x2) > max_level:
                continue

        if p_value < alpha:
            significance = [
                symbols[alpha_level_index]
                for alpha_level_index, alpha_level in enumerate(significance_levels)
                if p_value < alpha_level
            ][0]
            significance_str = "".join(significance)
            print(p_value, significance, significance_str)
        else:
            if print_ns:
                significance_str = "ns"
            else:
                significance_str = ""

        if len(significance_str) > 0:
            if shift_asterisks:
                y_val_mod = y_val + (i + 1) * y_range * 0.02
            else:
                y_val_mod = y_val + y_range * 0.02
            ax.text(
                x=(x1 + x2) / 2, y=y_val_mod, s=significance_str, ha="center", va="top"
            )
            if plot_hlines:
                line = Line2D(
                    [x1 + 0.1, x2 - 0.1], [y_val_mod, y_val_mod], lw=0.5, color="k"
                )
                ax.add_line(line)
                line.set_clip_on(False)
            # print((x1 + x2) / 2, y)


### Tools for plotting epigenetic profiles: 
from scipy.ndimage import gaussian_filter
def plot_profile(
    stack,
    background,
    idx_selected,
    ax=None,
    color="#000470",
    label="down",
    ylabel="Z-score",
    suptitle="",
    run_zscore=False,
    flank=200_000,
    resolution=10_000,
    ticks_step=50_000,
    abline_y=False,
    abline_x=False,
    scatter=True,
    return_source_data=False,
    run_smooth=True,
):
    """
    stack: numpy 2d array, full stack for the whole genome
    background: numpy 1d array, background values provided for normalization
    idx_selected: subset of the indexes of stack to plot & average
    ax: axes to put the plot to; if None, the whole new figure will be created
    """

    y_mean = np.nanmean(stack[idx_selected, :], axis=0)
    y_std = np.nanstd(y_mean)

    if run_smooth:
        y_smooth = gaussian_filter(y_mean, sigma=1.1)
        y_std_smooth = gaussian_filter(y_std, sigma=1.1)
    else:
        y_smooth = y_mean
        y_std_smooth = y_std

    x = np.arange(len(y_mean)) + 1

    if run_zscore:
        m, s = np.nanmean(background), np.nanstd(background)
        y_mean -= m
        y_mean /= s
        y_std /= s
        if run_smooth:
            y_smooth = gaussian_filter(y_mean, sigma=1.1)
            y_std_smooth = gaussian_filter(y_std, sigma=1.1)
        else:
            y_smooth = y_mean
            y_std_smooth = y_std

    if return_source_data:
        # y_max = np.nanmax( stack[idx_selected, :], axis=0)
        # y_min = np.nanmin( stack[idx_selected, :], axis=0)
        return x, y_mean, y_smooth, y_std_smooth, np.nanmean(background)

    if ax is None:
        f, axes = plt.subplots(1, 1, figsize=[5, 5])
        ax = axes

    if scatter:
        ax.scatter(x=x, y=y_mean, color=color, s=10)
    ax.plot(x, y_smooth, color=color, label=label)
    ax.fill_between(
        x, y_smooth - y_std_smooth, y_smooth + y_std_smooth, alpha=0.3, color=color
    )

    ax.set_xticks(np.arange(0, 2 * flank // resolution + 1, ticks_step // resolution))
    ax.set_xticklabels((np.arange(0, 2 * flank + 1, ticks_step) - flank) // 1_000)

    ax.set_xlim(0, 2 * flank // resolution)

    if abline_y:
        ax.plot(
            ax.get_xlim(),
            [np.nanmean(background), np.nanmean(background)],
            zorder=-1,
            linestyle="--",
            lw=1,
            color=color,
            alpha=0.8,
        )

    if abline_x:
        ax.plot(
            [0, 0],
            ax.get_ylim(),
            zorder=-1,
            linestyle="--",
            lw=1,
            color=color,
            alpha=0.8,
        )

    ax.set_title(suptitle)
    return ax


# Tools for plotting Hi-C pileups together:
def plot_pileups(
    stacks_list,
    idxs,
    titles,
    nrows=2,
    ncols=2,
    flank=200_000,
    resolution=10_000,
    figsize=(10, 10),
    vmin=-0.2,
    vmax=0.2,
    reduce_ticks=False,
    cmap="coolwarm",
    avscore=None,  # Otherwise the vector of dataframe (indexed as stacks) to aggregate
    add_counts=False,
    return_source_data=False,
    func=np.nanmean,
    f=None,
    axes=None,
):
    """Aggregate selected pileups, visualize matrices, and optionally return source data."""
    if not return_source_data:
        if axes is None and f is None:
            f, axes = plt.subplots(nrows, ncols, figsize=figsize)

    source_data = []
    for i, (idx_selected, title) in enumerate(zip(idxs, titles)):

        # Average pileup
        pile = func(stacks_list[list(idx_selected), :, :], axis=0)
        pile = fontanka.reflect(pile)
        score = func(avscore.loc[idx_selected])

        source_data.append([title, score, pile])

        if return_source_data:
            continue

        ax = axes.flatten()[i]

        im = ax.imshow(
            np.log2(pile), vmax=vmax, vmin=vmin, cmap=cmap, interpolation=None
        )

        if avscore is not None:
            ax.text(
                0.90,
                0.90,
                f"{score:.2f}",
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
                fontsize=20 if reduce_ticks else 17,
            )

        ticks_pixels = np.linspace(0, flank * 2 // resolution, 5)
        ax.set_xticks(ticks_pixels)
        ax.set_yticks(ticks_pixels)

        if not reduce_ticks:
            ticks_kbp = (
                (ticks_pixels - ticks_pixels[-1] / 2) * resolution // 1000
            ).astype(int)
            ax.set_xticklabels(ticks_kbp)
            ax.set_yticklabels(ticks_kbp)
            ax.set_xlabel("relative position, kbp")
            ax.set_ylabel("relative position, kbp")

        if reduce_ticks:
            ticks_kbp = ["" for x in ticks_pixels]
            ax.set_xticklabels(ticks_kbp)
            ax.set_yticklabels(ticks_kbp)

        if not reduce_ticks:
            divider = make_axes_locatable(ax)
            ax_colorbar = divider.append_axes("right", size="5%", pad=0.5)
            cbar = plt.colorbar(im, cax=ax_colorbar)

        if add_counts:
            title = f"{title}\n{len(idx_selected)}"
        ax.set_title(title)

    if return_source_data:
        return source_data

    f.tight_layout()

    return f


# Tools for plotting stratified genomic regions:
from sklearn.mixture import GaussianMixture
from scipy.stats import multivariate_normal


def fit_gaussian(data, xs):
    X = data.loc[:, xs].dropna().values
    idx = data.loc[:, xs].dropna().index.values

    # Fit 2D gaussian distribution with scikit learn
    model = GaussianMixture(
        n_components=1, covariance_type="full", max_iter=20, random_state=0
    )
    model.fit(X)

    return model.means_[0], model.covariances_[0]


def get_fdr(data, xs, mean, cov):
    X = data.loc[:, xs].dropna().values
    idx = data.loc[:, xs].dropna().index.values

    rv = multivariate_normal(mean=mean, cov=cov, allow_singular=False)

    pvs_collapsed = rv.pdf(X[:, :])

    pvs = pd.Series(np.zeros(len(data)) * np.nan)
    pvs.index = data.index
    pvs.loc[idx] = pvs_collapsed

    return pvs


def plot_2D_scatterplot(
    data,
    x,
    y,
    idxs,  # list of indexes for plotting
    titles=None,
    x_rng=(2, 6),
    y_rng=(2, 6),
    vmin_color=1,  # Color bars for counts
    vmax_color=50,
    xhist_stepsize=0.01,
    yhist_stepsize=0.01,
    colors=["blue"],
    cmaps=["blues"],
    ndots=120,
    figsize=(10, 7),
    xlabel="x",
    ylabel="y",
    cbar_label="counts",
    kind="dsshow",  # can be "scatter" or "regplot"
):

    norm = mpl.colors.PowerNorm(vmin=vmin_color, vmax=vmax_color, gamma=1 / 3)

    f, ax = plt.subplots(figsize=figsize)

    divider = make_axes_locatable(ax)
    ax_summary_x = divider.append_axes("bottom", size="20%", pad=0.1, sharex=ax)
    ax_summary_y = divider.append_axes("right", size="20%", pad=0.1, sharey=ax)
    density = False

    mask = (data[x] > -np.inf) & (data[y] > -np.inf)

    for i, idx in enumerate(idxs):
        if kind == "dsshow":
            im = dsshow(
                data.loc[idx, :],
                ds.Point(x, y),
                x_range=x_rng,
                y_range=y_rng,
                plot_height=ndots,
                plot_width=ndots,
                aspect="auto",
                cmap=cmaps[i],
                norm=norm,
                ax=ax,
            )
        elif kind == "scatter" or kind == "scatterplot":
            sns.scatterplot(
                data=data.loc[idx, :],
                x=x,
                y=y,
                color=colors[i],
                label=titles[i] if titles is not None else "",
                s=8,
                ax=ax,
            )
            ax.set_xlim(*x_rng)
            ax.set_ylim(*y_rng)
        elif (
            kind == "regplot" or kind == "reg"
        ):  # Plot regression line only for 'all' idxs

            if titles is not None:
                if titles[i] in ["all other", "all"]:

                    x_vals = data.loc[idx, x]
                    y_vals = data.loc[idx, y]
                    mask = np.isfinite(x_vals) & np.isfinite(y_vals)
                    r, p = scipy.stats.pearsonr(x_vals[mask], y_vals[mask])
                    print(f"PearsonR = {r:.3f}, pval = {p:.1e}; Lin Reg:")

                    model = sm.OLS(y_vals[mask], x_vals[mask]).fit()
                    print(model.summary())

                    sns.regplot(
                        data=data.loc[idx, :],
                        x=x,
                        y=y,
                        color=colors[i],
                        label=titles[i] if titles is not None else "",
                        scatter=False,
                        ax=ax,
                    )
                else:
                    sns.scatterplot(
                        data=data.loc[idx, :],
                        x=x,
                        y=y,
                        color=colors[i],
                        label=titles[i] if titles is not None else "",
                        s=8,
                        ax=ax,
                    )
            else:
                sns.scatterplot(
                    data=data.loc[idx, :],
                    x=x,
                    y=y,
                    color=colors[i],
                    label=titles[i] if titles is not None else "",
                    s=8,
                    ax=ax,
                )
            ax.set_xlim(*x_rng)
            ax.set_ylim(*y_rng)

        # Upper right hist:
        ax_summary_x.hist(
            data.loc[idx, :][x],
            bins=np.arange(*x_rng, xhist_stepsize),
            color=colors[i],
            histtype="step",
            density=density,
        )

        ax_summary_y.hist(
            data.loc[idx, :][y],
            bins=np.arange(*y_rng, yhist_stepsize),
            color=colors[i],
            histtype="step",
            density=density,
            orientation="horizontal",
        )
        if kind == "dsshow":
            ax_colorbar = divider.append_axes("right", size="5%", pad=0.1)
            cbar = plt.colorbar(
                im, cax=ax_colorbar, format=mpl.ticker.FuncFormatter(lambda x, pos: "")
            )
            ax_colorbar.set_ylabel(cbar_label)
        # else:
        #     ax.legend()

    # Major hist:
    ax_summary_x.hist(
        data[mask][x],
        bins=np.arange(*x_rng, xhist_stepsize),
        color="grey",
        histtype="step",
        density=density,
    )

    ax_summary_y.hist(
        data[mask][y],
        bins=np.arange(*y_rng, yhist_stepsize),
        color="grey",
        histtype="step",
        density=density,
        orientation="horizontal",
    )

    # Other
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax_summary_y.get_yticklabels(), visible=False)

    ax.set_ylabel(ylabel)
    ax_summary_x.set_xlabel(xlabel)

    return f, ax


def plot_2D_scatterplot_4way(
    data,
    x,
    y,
    x_rng=(2, 6),
    y_rng=(2, 6),
    vmin_color=1,  # Color bars for counts
    vmax_color=50,
    xhist_stepsize=0.01,
    yhist_stepsize=0.01,
    colors=[
        "blue",
        "purple",
        "black",
        "orange",
    ],  # listed from top left to bottom right np 2d flattening style
    cmaps=[
        "blues",
        "purples",
        "gray_r",
        "oranges",
    ],  # clockwise starting from upper left corner
    ndots=120,
    perc_x=70,
    perc_y=50,
    th_x=None,
    th_y=None,
    figsize=(10, 7),
    xlabel="x",
    ylabel="y",
    cbar_label="counts",
):

    norm = mpl.colors.PowerNorm(vmin=vmin_color, vmax=vmax_color, gamma=1 / 3)

    f, ax = plt.subplots(figsize=figsize)

    mask = (data[x] > -np.inf) & (data[y] > -np.inf)

    th_x = np.nanpercentile(data.loc[mask, x], perc_x)
    th_y = np.nanpercentile(data.loc[mask, y], perc_y)

    idx_top_both = data[mask & (data[x] > th_x) & (data[y] > th_y)].index
    idx_top_x = data[mask & (data[x] > th_x) & (data[y] <= th_y)].index
    idx_top_y = data[mask & (data[x] <= th_x) & (data[y] > th_y)].index
    idx_bottom = data[mask & (data[x] <= th_x) & (data[y] <= th_y)].index

    im_top_both = dsshow(
        data.loc[idx_top_both, :],
        ds.Point(x, y),
        x_range=x_rng,
        y_range=y_rng,
        plot_height=ndots,
        plot_width=ndots,
        aspect="auto",
        cmap=cmaps[1],
        norm=norm,
        ax=ax,
    )

    im_top_x = dsshow(
        data.loc[idx_top_x, :],
        ds.Point(x, y),
        x_range=x_rng,
        y_range=y_rng,
        plot_height=ndots,
        plot_width=ndots,
        aspect="auto",
        cmap=cmaps[3],
        norm=norm,
        ax=ax,
    )

    im = dsshow(
        data.loc[idx_bottom, :],
        ds.Point(x, y),
        x_range=x_rng,
        y_range=y_rng,
        plot_height=ndots,
        plot_width=ndots,
        aspect="auto",
        cmap=cmaps[2],
        norm=norm,
        ax=ax,
    )

    im_top_y = dsshow(
        data.loc[idx_top_y, :],
        ds.Point(x, y),
        x_range=x_rng,
        y_range=y_rng,
        plot_height=ndots,
        plot_width=ndots,
        aspect="auto",
        cmap=cmaps[0],
        norm=norm,
        ax=ax,
    )

    divider = make_axes_locatable(ax)
    ax_summary_x = divider.append_axes("bottom", size="20%", pad=0.1, sharex=ax)
    ax_summary_y = divider.append_axes("right", size="20%", pad=0.1, sharey=ax)

    density = False

    # Major hist:
    ax_summary_x.hist(
        data[mask][x],
        bins=np.arange(*x_rng, xhist_stepsize),
        color="grey",
        histtype="step",
        density=density,
    )

    ax_summary_y.hist(
        data[mask][y],
        bins=np.arange(*y_rng, yhist_stepsize),
        color="grey",
        histtype="step",
        density=density,
        orientation="horizontal",
    )

    # Upper right hist:
    ax_summary_x.hist(
        data.loc[idx_top_both, :][x],
        bins=np.arange(*x_rng, xhist_stepsize),
        color=colors[1],
        histtype="step",
        density=density,
    )

    ax_summary_y.hist(
        data.loc[idx_top_both, :][y],
        bins=np.arange(*y_rng, yhist_stepsize),
        color=colors[1],
        histtype="step",
        density=density,
        orientation="horizontal",
    )

    # Upper left:
    ax_summary_x.hist(
        data.loc[idx_top_x, :][x],
        bins=np.arange(*x_rng, xhist_stepsize),
        color=colors[3],
        histtype="step",
        density=density,
    )

    ax_summary_y.hist(
        data.loc[idx_top_x, :][y],
        bins=np.arange(*y_rng, yhist_stepsize),
        color=colors[3],
        histtype="step",
        density=density,
        orientation="horizontal",
    )

    # Lower left:
    ax_summary_x.hist(
        data.loc[idx_bottom, :][x],
        bins=np.arange(*x_rng, xhist_stepsize),
        color=colors[2],
        histtype="step",
        density=density,
    )

    ax_summary_y.hist(
        data.loc[idx_bottom, :][y],
        bins=np.arange(*y_rng, yhist_stepsize),
        color=colors[2],
        histtype="step",
        density=density,
        orientation="horizontal",
    )

    # Lower right:
    ax_summary_x.hist(
        data.loc[idx_top_y, :][x],
        bins=np.arange(*x_rng, xhist_stepsize),
        color=colors[0],
        histtype="step",
        density=density,
    )

    ax_summary_y.hist(
        data.loc[idx_top_y, :][y],
        bins=np.arange(*y_rng, yhist_stepsize),
        color=colors[0],
        histtype="step",
        density=density,
        orientation="horizontal",
    )

    # Other
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax_summary_y.get_yticklabels(), visible=False)

    ax_colorbar = divider.append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(
        im_top_both, cax=ax_colorbar, format=mpl.ticker.FuncFormatter(lambda x, pos: "")
    )

    ax_colorbar1 = divider.append_axes("right", size="5%", pad=0, sharey=ax_colorbar)
    plt.colorbar(
        im_top_x, cax=ax_colorbar1, format=mpl.ticker.FuncFormatter(lambda x, pos: "")
    )

    ax_colorbar2 = divider.append_axes("right", size="5%", pad=0, sharey=ax_colorbar)
    plt.colorbar(
        im, cax=ax_colorbar2, format=mpl.ticker.FuncFormatter(lambda x, pos: "")
    )

    ax_colorbar3 = divider.append_axes("right", size="5%", pad=0)
    plt.colorbar(im_top_y, cax=ax_colorbar3)

    ax.set_ylabel(ylabel)
    ax_summary_x.set_xlabel(xlabel)
    ax_colorbar3.set_ylabel(cbar_label)

    return f, (idx_top_both, idx_top_x, idx_bottom, idx_top_y, th_x, th_y)

### Tools for plotting Hi-C pileups:
def plot_heatmap_three_scales(pile, 
                 vmin_scale1=-0.5, 
                 vmax_scale1=0.5, 
                 vmin_scale2=-2.5, 
                 vmax_scale2=2.5, 
                 flank=200000, 
                 resolution=10000, 
                 interpolation="none", 
                 withlabels=True,
                 title=None,
                 cmap='RdBu_r',
                 cbar=True, 
                 log=False,
                 ax_hmap=None, 
                 grid_bar=None, 
                 fig=None):
    
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    if (ax_hmap is None) and (grid_bar is None) and (fig is None):
        print("Creating new figure...")
        fig = plt.figure()

        gs0 = gridspec.GridSpec(1, 2, width_ratios=[1, 0.1], height_ratios=[1], hspace=0)
        gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0], width_ratios=[1], height_ratios=[1])

        ax = plt.Subplot(fig, gs00[0])
        fig.add_subplot(ax)

        # grid spec for two right axes:
        gs01 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[1], height_ratios=[0.25, 1.0, 0.25], hspace=0)
        if cbar:
            ax_cbar1 = plt.Subplot(fig, gs01[0])
            ax_cbar2 = plt.Subplot(fig, gs01[1])
            ax_cbar3 = plt.Subplot(fig, gs01[2])
            fig.add_subplot(ax_cbar1)
            fig.add_subplot(ax_cbar2)
            fig.add_subplot(ax_cbar3)
        
    else: 
        if cbar:
            inner_grid = grid_bar.subgridspec(3, 1, height_ratios=[0.25, 1.0, 0.25], hspace=0)
            ax_cbar1, ax_cbar2, ax_cbar3 = inner_grid.subplots()
        
        ax = ax_hmap


    pile_1 = pile.copy()
    pile_1[(pile_1-vmax_scale1>=0) | (pile_1-vmin_scale1<0)] = np.nan

    pile_2 = pile.copy()
    pile_2[pile_2-vmax_scale1<0] = np.nan

    pile_3 = pile.copy()
    pile_3[pile_3-vmin_scale1>=0] = np.nan

    if log:
        norm_1 = mpl.colors.LogNorm(vmin_scale1, vmax_scale1)
    else:
        norm_1 = mpl.colors.Normalize(vmin_scale1, vmax_scale1)
    cmap_1 = cmap

    if log:
        norm_2 = mpl.colors.LogNorm(vmax_scale1, vmax_scale2)
    else:
        norm_2 = mpl.colors.Normalize(vmax_scale1, vmax_scale2)
    cmap_2 = mpl.colors.LinearSegmentedColormap.from_list("", [plt.cm.get_cmap(cmap, 1)(1), 'black'])

    if log:
        norm_3 = mpl.colors.LogNorm(vmin_scale2, vmin_scale1)
    else:
        norm_3 = mpl.colors.Normalize(vmin_scale2, vmin_scale1)
    cmap_3 = mpl.colors.LinearSegmentedColormap.from_list("", ['black', plt.cm.get_cmap(cmap, 2)(0)])

    # Visualize main map:
    im1 = ax.imshow(pile_1, cmap=cmap_1, norm=norm_1, interpolation=interpolation)
    im2 = ax.imshow(pile_2, cmap=cmap_2, norm=norm_2, interpolation=interpolation)
    im3 = ax.imshow(pile_3, cmap=cmap_3, norm=norm_3, interpolation=interpolation)

    if cbar:
        cbar1 = plt.colorbar(im1, cax=ax_cbar2)
        cbar2 = plt.colorbar(im2, cax=ax_cbar1)
        cbar3 = plt.colorbar(im3, cax=ax_cbar3)
        
        cbar1.ax.yaxis.set_ticks_position('left')
        cbar2.ax.yaxis.set_ticks_position('left')
        cbar3.ax.yaxis.set_ticks_position('left')
        
        cbar1.ax.spines[:].set_visible(False)
        cbar2.ax.spines[:].set_visible(False)
        cbar3.ax.spines[:].set_visible(False)

   
    ticks_pixels = np.linspace(0, flank*2//resolution,5)
    ticks_kbp = ((ticks_pixels-ticks_pixels[-1]/2)*resolution//1000).astype(int)
    if withlabels:
        ax.set_xlabel('relative position, kbp')
        ax.set_ylabel('relative position, kbp')
        
        ax.set_xticks(ticks_pixels)
        ax.set_yticks(ticks_pixels)
        ax.set_xticklabels(ticks_kbp)
        ax.set_yticklabels(ticks_kbp)
    else:
        ax.set_xticks(ticks_pixels)
        ax.set_yticks(ticks_pixels)
        ax.set_xticklabels([])
        ax.set_yticklabels([])        

    if title is not None:
        ax.set_title(title)
    # fig.tight_layout()
    return(fig)

def plot_heatmap(pile, 
                 vmin=-0.5, 
                 vmax_scale1=0.5, 
                 vmax_scale2=2.5, 
                 flank=200000, 
                 resolution=10000, 
                 interpolation="none", 
                 withlabels=True,
                 cmap='RdBu_r',
                 cbar=True, 
                 log=False,
                 ax_hmap=None, 
                 grid_bar=None, 
                 title=None,
                 fig=None):
    
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    if (ax_hmap is None) and (grid_bar is None) and (fig is None):
        print("Creating new figure...")
        fig = plt.figure()

        gs0 = gridspec.GridSpec(1, 2, width_ratios=[1, 0.1], height_ratios=[1], hspace=0)
        gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0], width_ratios=[1], height_ratios=[1])

        ax = plt.Subplot(fig, gs00[0])
        fig.add_subplot(ax)

        # grid spec for two right axes:
        gs01 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[1], height_ratios=[0.66, 0.33], hspace=0)
        if cbar:
            ax_cbar1 = plt.Subplot(fig, gs01[0])
            ax_cbar2 = plt.Subplot(fig, gs01[1])
            fig.add_subplot(ax_cbar1)
            fig.add_subplot(ax_cbar2)
        
    else: # Note that both ax_hmap and grid_bar should be supplied simultaneously

        if cbar:
            inner_grid = grid_bar.subgridspec(2, 1, height_ratios=[0.66, 0.33], hspace=0)
            ax_cbar1, ax_cbar2 = inner_grid.subplots()
        
        ax = ax_hmap


    pile_1 = pile.copy()
    pile_1[pile_1-vmax_scale1>=0] = np.nan

    pile_2 = pile.copy()
    pile_2[pile_2-vmax_scale1<0] = np.nan


    if log:
        norm_1 = mpl.colors.LogNorm(vmin, vmax_scale1)
    else:
        norm_1 = mpl.colors.Normalize(vmin, vmax_scale1)
    cmap_1 = cmap

    if log:
        norm_2 = mpl.colors.LogNorm(vmax_scale1, vmax_scale2)
    else:
        norm_2 = mpl.colors.Normalize(vmax_scale1, vmax_scale2)
    cmap_2 = mpl.colors.LinearSegmentedColormap.from_list("", [plt.cm.get_cmap(cmap, 1)(1), 'black'])

    # Visualize main map:
    im1 = ax.imshow(pile_1, cmap=cmap_1, norm=norm_1, interpolation=interpolation)
    im2 = ax.imshow(pile_2, cmap=cmap_2, norm=norm_2, interpolation=interpolation)

    if cbar:
        cbar1 = plt.colorbar(im1, cax=ax_cbar2)
        cbar2 = plt.colorbar(im2, cax=ax_cbar1)
        
        cbar1.ax.yaxis.set_ticks_position('left')
        cbar2.ax.yaxis.set_ticks_position('left')
        
        cbar1.ax.spines[:].set_visible(False)
        cbar2.ax.spines[:].set_visible(False)

   
    ticks_pixels = np.linspace(0, flank*2//resolution,5)
    ticks_kbp = ((ticks_pixels-ticks_pixels[-1]/2)*resolution//1000).astype(int)
    if withlabels:
        ax.set_xlabel('relative position, kbp')
        ax.set_ylabel('relative position, kbp')
        
        ax.set_xticks(ticks_pixels)
        ax.set_yticks(ticks_pixels)
        ax.set_xticklabels(ticks_kbp)
        ax.set_yticklabels(ticks_kbp)
    else:
        ax.set_xticks(ticks_pixels)
        ax.set_yticks(ticks_pixels)
        ax.set_xticklabels([])
        ax.set_yticklabels([])        

    if title is not None:
        ax.set_title(title)
        
    # fig.tight_layout()
    return(fig)
