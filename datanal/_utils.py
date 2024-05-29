"""Utility functions for the datanal package."""

from itertools import cycle
from typing import Any, Sequence
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import ticker
from matplotlib.axes import Axes


def auto_style(
    rc_mplstyle: dict[str, Any] | None = None,
    fname_mplstyle: str | None = None,
    palette_snsstyle: str | Sequence[str] | None = None,
    n_colors: int | None = None,
) -> None:
    """Initialize the matplotlib and seaborn styles.

    Parameters
    ----------
    rc_mplstyle : dict[str, Any] | None
        The matplotlib style for rcParams.
    fname_mplstyle : str | None
        The matplotlib style for figure.
    palette_snsstyle : str | Sequence[str] | None
        The seaborn palette style.
    n_colors : int | None
        The number of colors to generate.

    Returns
    -------
    None

    Notes
    -----
    Please put this function on the top of the script to enable global settings.
    """

    plt.style.use("datanal.GraphPadPrism")
    sns.set_palette("bright", n_colors=n_colors)

    if palette_snsstyle is not None:
        sns.set_palette(palette_snsstyle, n_colors=n_colors)
    if fname_mplstyle is not None:
        plt.style.use(fname_mplstyle)
    if rc_mplstyle is not None:
        plt.style.use(rc_mplstyle)


def auto_ticks(
    ax: Axes,
    *,
    left: float | None = None,
    right: float | None = None,
    bottom: float | None = None,
    top: float | None = None,
) -> None:
    """Set the major and minor ticks of an axis automatically.

    Parameters
    ----------
    ax : Axes
        The axis to be set.
    left : float | None, optional
        The left limit of the x-axis, by default None.
    right : float | None, optional
        The right limit of the x-axis, by default None.
    bottom : float | None, optional
        The bottom limit of the y-axis, by default None.
    top : float | None, optional
        The top limit of the y-axis, by default None.

    Returns
    -------
    None

    Notes
    -----
    Please put this function after the data are passed to the axis.
    """

    # Set the major and minor ticks of the x-axis
    if left is not None:
        ax.set_xlim(left=left)
    if right is not None:
        ax.set_xlim(right=right)

    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6, min_n_ticks=3))
    x_tick_num = len(ax.get_xticks())
    if x_tick_num >= 6:
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    elif x_tick_num >= 5:
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    else:
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

    # For viewing incomplete ticks
    if left is None and (ax.get_xticks()[0] != ax.get_xlim()[0]):
        ax.set_xlim(left=ax.get_xticks()[0])
    if right is None and (ax.get_xticks()[-1] != ax.get_xlim()[-1]):
        ax.set_xlim(right=ax.get_xticks()[-1])

    # Set the major and minor ticks of the y-axis
    if bottom is not None:
        ax.set_ylim(bottom=bottom)
    if top is not None:
        ax.set_ylim(top=top)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6, min_n_ticks=3))
    y_tick_num = len(ax.get_yticks())
    if y_tick_num >= 6:
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    elif y_tick_num >= 5:
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    else:
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))

    # For viewing incomplete ticks
    if bottom is None and (ax.get_yticks()[0] != ax.get_ylim()[0]):
        ax.set_ylim(bottom=ax.get_yticks()[0])
    if top is None and (ax.get_yticks()[-1] != ax.get_ylim()[-1]):
        ax.set_ylim(top=ax.get_yticks()[-1])


def auto_units(quant: float, old_unit: str = "M") -> tuple[str, float]:
    """Convert a quantity to a more appropriate unit.

    Parameters
    ----------
    quant : float
        The quantity to be converted.
    old_unit : str, optional
        The unit of the quantity, by default "M".

    Returns
    -------
    new_unit, factor : tuple[str, float]
        The converted unit and the conversion factor.
    """

    units = {
        "M": (1, 2e-1),
        "mM": (1e-3, 2e-4),
        "μM": (1e-6, 2e-7),
        "nM": (1e-9, 2e-10),
        "pM": (1e-12, 2e-13),
        "fM": (1e-15, 2e-16),
    }

    if old_unit not in units:
        warn(f"Unknown unit {old_unit} was detected.")
        return old_unit, 1

    old_scale = units[old_unit][0]
    quant *= old_scale

    return next(
        (
            (new_unit, old_scale / new_scale)
            for new_unit, (new_scale, bound) in units.items()
            if quant >= bound
        ),
        (old_unit, 1),
    )
