import os
import argparse
from warnings import warn

import numpy as np
import pandas as pd

# import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit

from ._utils import auto_style, auto_ticks, auto_units

class Enzyme:
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
    def sigmoidal(I: np.ndarray, IC50: float) -> np.ndarray:
        """The sigmoidal equation for enzyme inhibition.

        Inhibition rate ranges from 0 to 100. Hill coefficient is 1.

        Parameters
        ----------
        I : np.ndarray
            The inhibitor concentration.
        IC50 : float
            The half-maximal inhibitory concentration.

        Returns
        -------
        np.ndarray
            The inhibition rate.
        """

        return 100 / (1 + IC50 / I)

    @staticmethod
    def sigmoidal_4_param(
        I: np.ndarray, IC50: float, I_max: float, I_min: float
    ) -> np.ndarray:
        """The sigmoidal equation for enzyme inhibition with 4 parameters.

        Inhibition rate ranges from `I_min` to `I_max`. Hill coefficient is 1.

        Parameters
        ----------
        I : np.ndarray
            The inhibitor concentration.
        IC50 : float
            The half-maximal inhibitory concentration.
        I_max : float
            The maximum inhibition rate.
        I_min : float
            The minimum inhibition rate.

        Returns
        -------
        np.ndarray
            The inhibition rate.
        """

        return I_min + (I_max - I_min) / (1 + IC50 / I)

    @staticmethod
    def hill(I: np.ndarray, IC50: float, n: float) -> np.ndarray:
        """The Hill equation for enzyme inhibition.

        Inhibition rate ranges from 0 to 100, Hill coefficient is `n`.

        Parameters
        ----------
        I : np.ndarray
            The inhibitor concentration.
        IC50 : float
            The half-maximal inhibitory concentration.
        n : float
            The Hill coefficient.

        Returns
        -------
        np.ndarray
            The inhibition rate.
        """

        return 100 / (1 + (IC50 / I) ** n)

    @staticmethod
    def hill_4_param(
        I: np.ndarray, IC50: float, n: float, I_max: float, I_min: float
    ) -> np.ndarray:
        """The Hill equation for enzyme inhibition with 4 parameters.

        Inhibition rate ranges from `I_min` to `I_max`. Hill coefficient is `n`.

        Parameters
        ----------
        I : np.ndarray
            The inhibitor concentration.
        IC50 : float
            The half-maximal inhibitory concentration.
        n : float
            The Hill coefficient.
        I_max : float
            The maximum inhibition rate.
        I_min : float
            The minimum inhibition rate.

        Returns
        -------
        np.ndarray
            The inhibition rate.
        """

        return I_min + (I_max - I_min) / (1 + (IC50 / I) ** n)
