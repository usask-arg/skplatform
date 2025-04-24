from typing import Tuple, List, Optional
from numpy.typing import NDArray
from pathlib import Path
from datetime import datetime
from xarray import DataArray

import numpy as np
import xarray as xr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


def plot_files(filenames: Tuple[List[Path], List[str]], image_name: Tuple[str, Path],
                 time_step: int = 1, along_step: int =1, across_step: int = 1,
                 extent: Optional[List[float]] = [-180, 180, -90, 90]):
    """
    Plots the location of the satellite and the locations of the measurements in the
    files passed in. To view just the satellite track set 'along_step' and 'across_step' to 1.

    Parameters
    ----------
    filenames : Tuple[str, Path, List[Path], List[str]]
            A list of Paths or strings locating the data to plot.
    image_name : Tuple[str, Path]
            The name and location of the image to save the plot using.
    time_step : int, optional
            The step size to use when iterating through the measurements. A 'time_step' of
            4 would mean that every fourth observation would be included in the plot. The
            default is 1.
    along_step : int, optional
            The step size to take when iterating over the 'along' axis. A value of 5 would
            mean that every 5th angle in the along track FOV would be plotted. The default
            is 1.
    across_step : int, optional
            The step size to take when iterating over the 'across' axis. A value of 5 would
            mean that every 5th angle in the across track FOV would be plotted. the default
            is 1.
    extent : Optional[List[float]], optional
            The extend to plot, in the form of a list:
                [minimum_longitude, maximum_longitude, minimum_latitude, maximum_latitude]
            If set to None matplotlib will be allowed to size the plot. The default is the
            entire globe, [-180, 180, -90, 90].
    """

    results = []

    for file in filenames:
        results.append(xr.load_dataset(file))

    all = xr.concat(results, 'time')

    plot_dataset(all, image_name, time_step, along_step, across_step, extent)


def plot_dataset(ds: xr.Dataset, image_name: Tuple[str, Path],
                 time_step: int = 1, along_step: int =1, across_step: int = 1,
                 extent: Optional[List[float]] = None, axes = None):
    """
    Plots the data in a single dataset.

    Parameters
    ----------
    ds : xr.Dataset
            The dataset contianing the data to be plotted.
    image_name : Tuple[str, Path]
            The name and location of the image to save the plot using. This will be ignored
            if an axes aregument is provided.
    time_step : int, optional
            The step size to use when iterating through the measurements. A 'time_step' of
            4 would mean that every fourth observation would be included in the plot. The
            default is 1.
    along_step : int, optional
            The step size to take when iterating over the 'along' axis. A value of 5 would
            mean that every 5th angle in the along track FOV would be plotted. The default
            is 1.
    across_step : int, optional
            The step size to take when iterating over the 'across' axis. A value of 5 would
            mean that every 5th angle in the across track FOV would be plotted. The default
            is 1.
    extent : Optional[List[float]], optional
            The extend to plot, in the form of a list:
                [minimum_longitude, maximum_longitude, minimum_latitude, maximum_latitude]
            If set to None matplotlib will be allowed to size the plot. The default is None.
    axes : _type_, optional
            If None a new figure and image is created, and the plot is made using them. If
            not None then an exiting axes is plotted to. The image will not be saved if
            this in not None. The default is None.

    Returns
    -------
    _type_ ax
        _description_
            Returns the axes that was either created or plotted to.
    """
    if axes is None:
        fig, ax = plt.subplots(1,1, figsize=(16, 9), dpi=150, subplot_kw={'projection':ccrs.PlateCarree()})
    else:
        ax = axes

    ax.coastlines()
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='-')
    if extent is not None:
        ax.set_extent(extent)

    ax.plot(ds.observer_longitude[::time_step], ds.observer_latitude[::time_step], 'ok')

    lat = ds.latitude[::time_step, ::along_step, ::across_step].to_numpy().flatten()
    lon = ds.longitude[::time_step, ::along_step, ::across_step].to_numpy().flatten()

    ax.scatter(lon, lat, s=0.1)

    # if axes is set ignore the image name (fig wasn't created anyway)
    if (image_name is not None) and (axes is None):
        fig.savefig(image_name)

    return ax


def plot_observation(ds: xr.Dataset, image_name: Optional[Tuple[str, Path]],
                 time: Tuple[List[int], List[str], int, str] = 0, along_step: int =1, across_step: int = 1,
                 extent: Optional[List[float]] = None, axes = None):
    """
    Plots a single 'observation' from a dataset, where an observation is all of the measurements
    taken at a single time.

    Parameters
    ----------
    ds : xr.Dataset
            The dataset contianing the observation to be plotted.
    image_name : Optional[Tuple[str, Path]]
            The name of the
    time : Tuple[List[int], List[str], int, str], optional
            The time of the observation to be plotted, either as a sting or as an integer index
            to the time array in the dataset. If a string, the closest dataset time to the passed
            in value will be selected. Single values can be used. The default is 0, the first
            time in the dataset.
    along_step : int, optional
            The step size to take when iterating over the 'along' axis. A value of 5 would
            mean that every 5th angle in the along track FOV would be plotted. The default is 1.
    across_step : int, optional
            The step size to take when iterating over the 'across' axis. A value of 5 would
            mean that every 5th angle in the across track FOV would be plotted. The default is 1.
    extent : Optional[List[float]], optional
            The extend to plot, in the form of a list:
                [minimum_longitude, maximum_longitude, minimum_latitude, maximum_latitude]
            If set to None matplotlib will be allowed to size the plot. The default is None.
    axes : _type_, optional
            If None a new figure and image is created, and the plot is made using them. If not None
            then an exiting axes is plotted to. The image will not be saved if this in not None. The
            default is None.

    Returns
    -------
    ax : matplotlib.axes.Axes
        Either the axis object created by this function or the axes that was passed in as an argument.
    """

    if axes is None:
        fig, ax = plt.subplots(1,1, figsize=(16, 9), dpi=150, subplot_kw={'projection':ccrs.PlateCarree()})
    else:
        ax = axes

    ax.coastlines()
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='-')
    if extent is not None:
        ax.set_extent(extent)

    if not isinstance(time, list):
        time = [time]

    for t in time:
        if isinstance(t, int):
            # an direct index
            selected = ds.isel(time=t, drop=True)
        elif isinstance(t, str):
            selected = ds.sel(time=t, drop=True, method='nearest')

        ax.plot(selected.observer_longitude, selected.observer_latitude, 'ok')

        lat = selected.latitude[::along_step, ::across_step].to_numpy().flatten()
        lon = selected.longitude[::along_step, ::across_step].to_numpy().flatten()

        ax.scatter(lon, lat, s=0.1)

    # if axes is set ignore the image name (fig wasn't created anyway)
    if (image_name is not None) and (axes is None):
        fig.savefig(image_name)

    return ax


def blank_plot():
    """
    Returns figure and axes objects setup for geographic plotting. The projection is set to PlateCarree,
    coast lines and grid lines are also configured.
    """
    fig, ax = plt.subplots(1,1, figsize=(16, 9), dpi=150, subplot_kw={'projection':ccrs.PlateCarree()})
    ax.coastlines()
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='-')
    return fig, ax


def convert_dataset_time_to_list(dataset_time: DataArray[np.datetime64]):
    """
    Converts time from a dataset (a DataArray of np.datetime64) into a list
    of time strings for use in selecting observations to plot. Useful when
    trying to select the same time in different files.

    Parameters
    ----------
    dataset_time : DataArray[np.datetime64]
        The time to convert to strings.

    Return
    ------
    List[str]
        A list of strings representing the times in thegit push `dataset_time` arguement.
    """
    return [str(t) for t in np.datetime_as_string(dataset_time.values)]

# can be called from the repo with one line
# python -c "from plotlib import plot_files_2; from pathlib import Path; plot_files_2([ file for file in Path('./output/orbits').glob('orbit_*')], image_name='one_day2.png' )"

def plot_dataset_satellite_path(ds: xr.Dataset, image_name: Tuple[str, Path],
                 time_step: int = 1, extent: Optional[List[float]] = None, axes = None):
    """
    Plots the only the satellite path from a single dataset. Measurement locations are not plotted.
    Each orbit will be plotted a different color, using the current default color cycle.

    Parameters
    ----------
    ds : xr.Dataset
            The dataset contianing the data to be plotted.
    image_name : Tuple[str, Path]
            The name and location of the image to save the plot using. This will be ignored
            if an axes aregument is provided.
    time_step : int, optional
            The step size to use when iterating through the measurements. A 'time_step' of
            4 would mean that every fourth observation would be included in the plot. The
            default is 1.
    extent : Optional[List[float]], optional
            The extend to plot, in the form of a list:
                [minimum_longitude, maximum_longitude, minimum_latitude, maximum_latitude]
            If set to None matplotlib will be allowed to size the plot. The default is None.
    axes : matplotlib.axes.Axes, optional
            If None a new figure and image is created, and the plot is made using them. If not None
            then an exiting axes is plotted to. The image will not be saved if this in not None. The
            default is None.

    Returns
    -------
    ax : matplotlib.axes.Axes
        Returns the axes that was either created or plotted to.
    """

    if axes is None:
        fig, ax = plt.subplots(1,1, figsize=(16, 9), dpi=150, subplot_kw={'projection':ccrs.PlateCarree()})
    else:
        ax = axes

    ax.coastlines()
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='-')
    if extent is not None:
        ax.set_extent(extent)

    for orbit in np.unique(ds.orbit.values):
        subset = ds.where(ds.orbit==orbit, drop=True)
        ax.plot(subset.observer_longitude[::time_step], subset.observer_latitude[::time_step])

    # if axes is set ignore the image name (fig wasn't created anyway)
    if (image_name is not None) and (axes is None):
        fig.savefig(image_name)

    return ax