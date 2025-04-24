from typing import Tuple

import logging
import numpy as np
import xarray as xr
from pathlib import Path

def convert_to_single_orbit_files(search_path: Tuple[Path, str], output_path: Tuple[Path, str],
                                  glob_string: str = '*', output_prefix: str = 'orbit'):
    """
    Finds all the files in the given directory that match the given pattern and converts them
    into files containing a single orbit each. The output files will be saved with the format
    output_prefix + '_###.nc'

    Parameters
    ----------
    search_path : Tuple[Path, str]
        Either a string or a Path object indicating the folder to search for files.
    output_path : Tuple[Path, str]
        Either a string or a Path object indicating the folder to save the output files in.
    glob_string : str, optional
        The pattern filenames must match to be included. The default, '*', matches all files.
    output_prefix: str, optional
        The string prepended to the name of the converted files. Th default is 'orbit'.
    """
    # this relies heavily on the filename and folder structure, but opening all the files would be a big
    # performance hit and possibly take too much memory if not lazy loading.
    # it also produces nonsense if output of two runs are mixed
    existing_ds = None
    if isinstance(search_path, str):
        input_dir = Path(search_path)
    elif isinstance(search_path, Path):
        input_dir = search_path
    else:
        logging.error('Could not convert search_path into a Path object.')
        return

    if isinstance(output_path, str):
        output_dir = Path(output_path)
    elif isinstance(output_path, Path):
        output_dir = output_path
    else:
        logging.error('Could not convert output_path into a Path object.')
        return

    files = input_dir.glob(glob_string)
    files = [file for file in files]
    files.sort()  #

    existing_ds = None
    for file in files:
        try:
            new_ds = xr.load_dataset(file)
            if existing_ds is not None:
                existing_ds = xr.concat([existing_ds, new_ds], 'time')
            else:
                existing_ds = new_ds
            # if the files are loaded in time order then everytime we have more than
            # one orbit number we should have all the data for the lowest orbit number loaded
            orbit_nums = np.unique(existing_ds.orbit)
            while orbit_nums.size > 1:
                netcdf_idx = np.where( np.isclose(existing_ds.orbit, orbit_nums[0]) )[0]
                lowest_orbit = existing_ds.isel(time=netcdf_idx)
                # update the time, otherwise it will use the time of the first file
                # for all of them
                lowest_orbit.attrs['start_time'] = lowest_orbit.time[0].dt.strftime('%Y-%m-%d %H:%M:%S').item()
                lowest_orbit.attrs['end_time'] = lowest_orbit.time[-1].dt.strftime('%Y-%m-%d %H:%M:%S').item()
                lowest_orbit.to_netcdf(output_dir.joinpath(f'{output_prefix}_{orbit_nums[0]:04}.nc'))
                keep_idx = np.where( ~np.isclose(existing_ds.orbit, orbit_nums[0]) )[0]
                existing_ds = existing_ds.isel(time=keep_idx, drop=True)
                orbit_nums = np.unique(existing_ds.orbit)
        except ValueError as ex:
            # this occurs when there is a non-netcdf file in the directory, ex. 'orbit_1.png'
            pass

    # export any remaining data
    if existing_ds is not None:
        existing_ds.attrs['start_time'] = existing_ds.time[0].dt.strftime('%Y-%m-%d %H:%M:%S').item()
        existing_ds.attrs['end_time'] = existing_ds.time[-1].dt.strftime('%Y-%m-%d %H:%M:%S').item()
        existing_ds.to_netcdf(output_dir.joinpath(f'{output_prefix}_{existing_ds.orbit[0].item():04}.nc'))

def convert_to_single_file(search_path: Tuple[Path, str], output_file: Tuple[Path, str],
                           glob_string: str = '*' ):
    """
    Finds all the files in the given directory that match the given pattern and converts them
    into a single file.

    Parameters
    ----------
    search_path : Tuple[Path, str]
        Either a string or a Path object indicating the folder to search for files.
    output_file : Tuple[Path, str]
        Either a string or a Path object indicating the folder and name of the output file.
    glob_string : str, optional
        The pattern filenames must match to be included. The default, '*', matches all files.
    """

    existing_ds = None
    if isinstance(search_path, str):
        input_dir = Path(search_path)
    elif isinstance(search_path, Path):
        input_dir = search_path
    else:
        logging.error('Could not convert search_path into a Path object.')
        return

    if isinstance(output_file, str):
        output_filepath = Path(output_file)
    elif isinstance(output_file, Path):
        output_filepath = output_file
    else:
        logging.error('Could not convert output_path into a Path object.')
        return

    files = input_dir.glob(glob_string)
    files = [file for file in files]
    files.sort()  #

    existing_ds = None
    for file in files:
        try:
            new_ds = xr.load_dataset(file)
            if existing_ds is not None:
                existing_ds = xr.concat([existing_ds, new_ds], 'time')
            else:
                existing_ds = new_ds
        except ValueError as ex:
            # this occurs when there is a non-netcdf file in the directory, ex. 'orbit_1.png'
            pass

    # export any remaining data
    if existing_ds is not None:
        existing_ds.attrs['start_time'] = existing_ds.time[0].dt.strftime('%Y-%m-%d %H:%M:%S').item()
        existing_ds.attrs['end_time'] = existing_ds.time[-1].dt.strftime('%Y-%m-%d %H:%M:%S').item()
        existing_ds.to_netcdf(output_filepath)