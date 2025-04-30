from typing import Tuple, List

import xarray as xr
import simplekml

from pathlib import Path

def satellite_path_to_kml(ds: xr.Dataset, save_as: Tuple[Path, str], time_step: int = 1, color: str = 'green'):
    """
    Generates a kml file for Google Earth that contains pins for the path of the satellite in the passed in dataset.
    To view multiple datasets together, call this function on multiple files and load all of the kml files. Google
    Earth can struggle when a large number of pins are loaded.

    Parameters
    ----------
    ds : xr.Dataset
            The dataset to plot.
    save_as : Tuple[Path, str]
            The name and path so save the kml as.
    time_step : int, optional
            The step size to use when iterating through the measurements. A value of 3, for instance, means
            every third measurement will be plotted. Useful for limiting the number of points drawn. The
            default value is 1.
    color : str, optional
            The color to use for the pins. Can be either a named color or a RGBA hex string, 'ffffff00' for
            example. The default color is green.
    """
    kml = simplekml.Kml()

    # the points accept hex values 'fffffffff' (RGBA)
    # look up the named colors and use it if found, else assume the color is a hex
    colors = [c for c in dir(simplekml.Color) if isinstance(getattr(simplekml.Color, c), str)]
    if color in colors:
        color_val = getattr(simplekml.Color, color)
    else:
        color_val = color

    for t in range(0, ds.time.shape[0], time_step):
        coords = (ds.observer_longitude[t].item(), ds.observer_latitude[t].item())
        point = kml.newpoint(name='', coords=[coords])
        point.style.iconstyle.color = color_val

    kml.save(save_as)

def satellite_observation_to_kml(ds, save_as: Tuple[Path, str], time: List[Tuple[int, str]] = 1, color: str = 'green'):
    """
    Generates a kml file for Google Earth that contains pins for the measurments of the satellite in the passed in dataset.
    To view multiple datasets together, call this function on multiple files and load all of the kml files. Google
    Earth can struggle when a large number of pins are loaded.

    Parameters
    ----------
    ds : _type_
            The dataset to plot.
    save_as : Tuple[Path, str]
            The name and path so save the kml as.
    time : List[Tuple[int, str]], optional
            Either a integer index or a time string (ex. '2019-08-01 00:10:04') used for selecting the measurements to draw.
            Can also be a list of integer indices or time strings. The list cannot be a mix of integers and strings. The
            default is 1.
    along_step : int, optional
            The step size to take when iterating over the 'along' axis of the measurements. A value of 5 would mean that
            every 5th along angle would be plotted. The default is 1.
    across_step : int, optional
            The step size to take when iterating over the 'across' axis of the measurements. A value of 5 would mean that
            every 5th across angle would be plotted. The default is 1.
    color : str, optional
            The color to use for the pins. Can be either a named color or a RGBA hex string, 'ffffff00' for
            example. The default is green.
    """
    kml = simplekml.Kml()
    # the points accept hex values 'fffffffff' (RGBA)
    # look up the named colors and use it if found, else assume the color is a hex
    colors = [c for c in dir(simplekml.Color) if isinstance(getattr(simplekml.Color, c), str)]
    if color in colors:
        color_val = getattr(simplekml.Color, color)
    else:
        color_val = color

    if not isinstance(time, list):
        time = [time]

    for t in time:
        if isinstance(t, int):
            # an direct index
            selected = ds.isel(time=t, drop=True)
        elif isinstance(t, str):
            selected = ds.sel(time=t, drop=True, method='nearest')

        for along in ds.along:
            for across in ds.across:
                coords = (selected.longitude[along, across].item(), selected.latitude[along, across].item())
                point = kml.newpoint(name='', coords=[coords])
                point.style.iconstyle.color = color_val

    kml.save(save_as)

