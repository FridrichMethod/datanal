import argparse
import os
from typing import Any
from warnings import warn

# import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import curve_fit

from datanal._utils import auto_style, auto_ticks, auto_units


class Enz:
    @staticmethod
    def michaelis_menten(S: np.ndarray, V_m: float, K_M: float) -> np.ndarray:
        """The Michaelis-Menten equation for enzyme kinetics.

        Parameters
        ----------
        S : np.ndarray
            The substrate concentration.
        V_m : float
            The maximum rate of the reaction.
        K_M : float
            The Michaelis constant.

        Returns
        -------
        np.ndarray
            The reaction rate.
        """

        return V_m * S / (K_M + S)

    @staticmethod
    def sigmoidal(
        I: np.ndarray,
        IC50: float,
        *,
        I_max: float = 100,
        I_min: float = 0,
        hill_coeff: float = 1,
        monotonic: int = 1,
    ) -> np.ndarray:
        """The Hill equation for enzyme inhibition with 4 parameters.

        Inhibition rate ranges from `I_min` to `I_max`. Hill coefficient is `hill_coeff`.

        Parameters
        ----------
        I : np.ndarray
            The inhibitor concentration.
        IC50 : float
            The half-maximal inhibitory concentration.
        I_max : float, optional
            The maximum inhibition rate, by default 100.
        I_min : float, optional
            The minimum inhibition rate, by default 0.
        hill_coeff : float, optional
            The Hill coefficient, by default 1.
        monotonic : int, optional
            The monotonicity of the curve, by default 1.

        Returns
        -------
        np.ndarray
            The inhibition rate.
        """

        if monotonic not in {-1, 1}:
            raise ValueError("Monotonicity must be either -1 or 1.")

        return I_min + (I_max - I_min) / (1 + (IC50 / I) ** (hill_coeff * monotonic))

    @staticmethod
    def michaelis_menten_fit(
        x_data: np.ndarray, y_data: np.ndarray
    ) -> tuple[float, float] | tuple[None, None]:
        """Fit the Michaelis-Menten equation to the data.

        Parameters
        ----------
        x_data : np.ndarray
            The substrate concentration.
        y_data : np.ndarray
            The reaction rate.

        Returns
        -------
        tuple[float, float] | tuple[None, None]
            The fitted parameters (V_m, K_M) or None if the fitting fails.
        """

        # Initial guess
        initial_guess = ...

        # Bounds
        lower_bounds = ...
        upper_bounds = ...

        # Fit affinity
        try:
            popt, pcov, *_ = curve_fit(
                Enz.michaelis_menten,
                x_data,
                y_data,
                p0=initial_guess,
                bounds=(lower_bounds, upper_bounds),
            )
        except RuntimeError:
            # Some affinity data are not fitted successfully
            return None, None

        # Check the standard error of the fitted parameters
        ...

        return popt, pcov

    @staticmethod
    def sigmoidal_fit(
        x_data: np.ndarray, y_data: np.ndarray
    ) -> tuple[float, float] | tuple[None, None]:
        """Fit the sigmoidal equation to the data.

        Parameters
        ----------
        x_data : np.ndarray
            The inhibitor concentration.
        y_data : np.ndarray
            The inhibition rate.

        Returns
        -------
        ...
        """

        # Initial guess
        initial_guess = ...

        # Bounds
        lower_bounds = ...
        upper_bounds = ...

        # Fit inhibition
        try:
            popt, pcov, *_ = curve_fit(
                Enz.sigmoidal,
                x_data,
                y_data,
                p0=initial_guess,
                bounds=(lower_bounds, upper_bounds),
            )
        except RuntimeError:
            # Some inhibition data are not fitted successfully
            return None, None

        # Check the standard error of the fitted parameters
        ...

        return popt, pcov


class EnzPlot:
    def __init__(
        self,
        affinity_file: str,
        *,
        rc_mplstyle: dict[str, Any] | None = None,
        fname_mplstyle: str | None = None,
        palette_snsstyle: str | list[str] | None = None,
    ) -> None:
        """Initialize the matplotlib and seaborn styles and the affinity file.

        Parameters
        ----------
        affinity_file : str
            The enzyme data file.
        rc_mplstyle : dict[str, Any] | None
            The matplotlib style for rcParams.
        fname_mplstyle : str | None
            The matplotlib style for figure.
        palette_snsstyle : str | list[str] | None
            The seaborn style for palette.

        Returns
        -------
        None
        """

        auto_style(rc_mplstyle, fname_mplstyle, palette_snsstyle)

        if os.path.exists(affinity_file):
            self._affinity_file = affinity_file
        else:
            raise FileNotFoundError(f"{affinity_file} does not exist.")

    def __repr__(self) -> str:
        """Return the name of the affinity file."""

        return os.path.basename(self._affinity_file)

    def __len__(self) -> int:
        """Return the number of data points in the affinity file."""

        return len(self.affinity_data)

    @property
    def affinity_data(self) -> pd.DataFrame:
        """Read the first two data frames (concentration and response) of the affinity file."""

        df = pd.read_excel(self._affinity_file)

        return ...
