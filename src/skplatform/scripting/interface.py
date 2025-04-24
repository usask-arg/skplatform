from typing import List, Tuple, Optional
from datetime import datetime, timedelta, timezone
import numpy as np
import xarray as xr
from xarray import Dataset
from pathlib import Path
import logging
from copy import deepcopy
from multiprocessing import Pool
from .instruments import AOSSky, TICFIRE, PropagatorInstrument


def read_time_string(time_str: str):
    """
    Converts a time string in the format 'YYYY-mm-dd HH:MM:SS', and returns
    its as a datetime object with a UTC timezone.
    """
    # the ValueError is message is pretty good
    time = datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')
    return time.replace(tzinfo=timezone.utc)


def default_tle():
    """
    The default TLE as a list of two strings. This will be updated as
    newer TLE's are avaiable. Currently taken from
    AOS_Sky_Descending_v2025-01-29_Epoch2019.tle. Provided so that
    scripts and tests have a default available that does not require
    locating and/or downloading a TLE file.
    """
    # taken from AOS_Sky_Descending_v2025-01-29_Epoch2019.tle
    # directly entered so don't rely on the file
    return ['1 99999U 99999 A  19213.00000000  .00000000  00000-0  00000-0 0    10\n',
            '2 99999 097.2222 331.8182 0001962 091.0782 359.0176 15.37843164    16\n']

def open_tle(tle: Path | str):
    """
    Open a TLE file and return it as a two-element list of strings.

    Parameters
    ----------
    tle: Path | str
        Either a string or a path indicating the file to read. The file
        must be valid. Files with a trailing newline will be read, but
        all other variations will not.

    Returns
    -------
    List[str]
        A 2 element list with the first line as the first element and the
        second line as the second element.
    """
    try:
        with open(tle) as file:
            tle_lines = file.readlines()
            # sometimes an extra newline is at the end of the file
            if len(tle_lines) > 2 and tle_lines[2] == '\n':
                tle_lines = tle_lines[:2]
    except IOError as ex:
        logging.error(f'There was an error with the TLE file {tle}.')
        raise
    return tle_lines

def simulate_orbit(
        instrument: str | PropagatorInstrument = 'aos-sky',
        output: Path | str = None,
        tle: Path | str = '',
        start: Optional[datetime | str] = None,
        end: Optional[datetime | str] = None,
        length: float = 60,
        report_epoch: bool = False,
        threads: int = 1) -> Optional[Dataset]:
    """
    The main interface function for propagating orbits.

    Parameters
    ----------
    instrument : str | PropagatorInstrument, optional
        by default 'aos-sky'
            The instrument to run. Can be a string identifiying the instrument or an instance of a subclass
            of PropagatorInstrument.

            If a string is used it must be either "AOS-Sky" or "TICFIRE". The names
            are not case sensitive.

            It is expected that any PropagatorInstrument instance passed in is fully
            initialized. The tle argument will not be used in this case.

            The default is 'aos-sky'.
    output : Optional[Path | str], optional
            The path and filename to save the simulation output with. Can take on 3
            different values:\n\n
            None: no file will be saved.\n\n
            A directory: the file will be saved with an name using the format:
            ``{instrument}_{starttime}_to_{endtime}.nc``\n\n
            An empty string: a file will be saved to the directory the script is running in, using automatic naming.\n\n
            A full name and path: a file will be saved with give name and path.\n\n
            Any existing file will be overwritten.
            The default is None.
    tle : Path | str, optional
            The file and path of a TLE to load. If not supplied the default TLE is used.
            The default is ''.
    start : Optional[datetime | str], optional
            The time to start the simulation from. This can be before the TLE epoch.
            The format is strict, and must be "YYYY-MM-DD hh:mm:ss", quotes included.
            If not set the simulation starts from the TLE epoch. The default is None.
    end : Optional[datetime | str], optional
            The time to end the simulation at. This can be before the TLE epoch, but
            only if --start is also provided. --start must be before --end. The
            format is strict, and must be "YYYY-MM-DD hh:mm:ss", quotes included.
            Overrides the --length argument. The default is None.
    length : float, optional
            The number of seconds to run the simulation for. Only used if --end is not set.
            The default is 60 seconds.
    report_epoch : bool, optional
            If set to 1, then the TLE epoch will be printed to the console and the
            program will exit. The default is False.
    threads : int, optional, by default 1
            The number of threads to use. If greater than 1 then the simulations changes to
            a multi-threaded mode (technically multi-process) from the default single
            threaded mode. When in multi-threading mode, the function must be called
            from inside an "if __name__ == '__main__':" statement, or the multithreading will
            fail. "output" must be a directory when in multithreading mode, and all output
            files will be saved using the default naming pattern.

            The time to simulate will be broken into chunks of 1024 measurements, ie. each
            thread will be given 1024*instrument frequency seconds to simulate. If there are
            more chunks than threads this argument sets the maximum number of threads that can
            exist at a time. Each time a thread completes it will save its results as an
            individual file. The 'conversion' module contains functions to process the output
            of this mode into a single file, or a single file per orbit.

            As a rule of thumb the number of threads should not exceed the number of physical
            cores available. If the computer is still being used for other tasks it is
            recommended to set threads to number of physical cores - 2 so that there are
            cores available for other tasks. Note that operating systems often report
            the number of logical cores, where each physical core with hyper-threading
            enabled is counted twice.

    Returns
    -------
    Optional[Dataset]
        Returns the dataset generated if ran in single thread mode. In multi-thread mode
        None is returned.
    """
    # load the TLE
    tle_lines = default_tle()
    if tle != '':
        tle_lines = open_tle(tle)

    instrument_name = _get_instrument_name(instrument)

    inst = _create_instrument(instrument, instrument_name, tle_lines)

    # setup time
    tle_time = inst.platform.platform_locator.time

    if report_epoch:
        print(f'The TLE epoch is {tle_time}.')
        return

    _set_instrument_time(inst, start, end, length, tle_time)

    ds = None

    if threads == 1:
        inst.propagate()
        ds = inst.get_dataset()
        if output is not None:
            auto_name = _make_filename_from_ds_time(instrument_name, ds)

            if output == '':
                ds.to_netcdf(auto_name)
            elif (isinstance(output, Path) and output.is_dir()):
                ds.to_netcdf(output.joinpath(auto_name))
            else:
                ds.to_netcdf(output)
    else:
        if output is None:
            logging.error('Argument output is None and threads is > 1. Either provide a directory '
                          'to save to in output or switch to single threaded mode (threads=1). Use an '
                          'empty string to save to the directory the script is being ran from.')
            raise TypeError

        if isinstance(output, str):
            output = Path(output)

        if output.is_file():
            logging.error('Trying to call a multi-threaded simulation with output as a file. '
                          'Mutlithreaded simulations only support output as a directory.')
            raise TypeError

        # split the simulation time into chunks and pass them off to a function
        # for processing. Each thread savves it's results to disk, avoiding
        # both trying to hold the entire orbit in memory and problems with
        # resource sharing until IO is overloaded.
        if isinstance(inst.frequency, float):
            chunk_size = 1024*inst.frequency
        elif isinstance(inst.frequency, (List, np.ndarray)):
            chunk_size = 1024 * np.sum(inst.frequency)

        start_times = np.arange(inst.start_time, inst.end_time, chunk_size)
        lengths = np.diff(start_times)
        # append last interval
        lengths = np.append(lengths, inst.end_time - start_times[-1])

        pool_args = [(deepcopy(instrument), deepcopy(tle_lines), start_times[k], lengths[k], output) for k in range(lengths.size)]
        # At least on Windows this works to prevent inifinte loops of crashing processes
        try:
            # use min so don't start too many threads, if need fewer than the threads argument
            with Pool(min(len(lengths), threads)) as pool:
                pool.starmap(_run_thread, pool_args)
        except RuntimeError as ex:
            logging.exception('A RuntimeError occured starting simulate_orbit() in multi-threaded mode.\n'
                              'Is simulate_orbit() protected in a "if __name__ == "__main__":" clause? '
                              'This is a requrement for running in multi-threaded mode.\n\n')
            #raise ex  # cannot re-raise, or it will not terminate
            return

    return ds

def _run_thread(instrument: str, tle_lines: List[str], start_seconds: float, length_seconds: float, output: str | Path):
    """
    Calculates the orbit over a portion of the entire simulation and saves
    the results to disk. Only intended to be used with `simulate_orbit()`.
    Inputs are assumed to have been validated by `simulate_orbit()`.

    Parameters
    ----------
    instrument : str | PropagatorInstrument
        A copy of a configured instrument or a string designateing either AOS-Sky or TICFIRE.
    tle_lines : List[str]
        The TLE of the orbit being simulated, as a list of two strings. Each string should be
        one line of the TLE.
    start_seconds : float
        The number of seconds offset from the TLE epoch to begin this portion of the
        simulation at. Can be negative.
    length_seconds : float
        The number of seconds to simulate the orbit for.
    output : str | Path
        The location to save the results at.
    """
    # reconstruct instrument inside the thread to avoid any problems with
    # shared objects. Strings and floats are immutable; instrument is only read from
    instrument_name = _get_instrument_name(instrument)

    inst = _create_instrument(instrument, instrument_name, tle_lines)

    # these should be calculated in the calling function
    inst.set_start_time(start_seconds)
    inst.set_end_time(length_seconds + start_seconds)

    inst.propagate()
    ds = inst.get_dataset()

    filename = _make_filename_from_ds_time(instrument_name, ds)

    output_path = Path(output)
    ds.to_netcdf(output_path.joinpath(filename))


def _get_instrument_name(instrument: str | PropagatorInstrument) -> str:
    """
    Either formats a string instrument name or tries to get the name of
    an instrument class instance.

    Parameters
    ----------
    instrument : str | PropagatorInstrument
        Either a string or an instance of the PropagatorInstrument or its
        subclasses. A string will be changed to lower case, while a class
        instance will be queried for it's class name. If the class has a
        "name" attribute it will be used, otherwise the name of the class
        in the code will be formatted and used.

    Returns
    -------
    str
        The name of the instrument, or the name of it's class.
    """
    if isinstance(instrument, str):
        instrument_name = instrument.lower()
    elif isinstance(instrument, PropagatorInstrument):
        if hasattr(instrument, 'name'):
            instrument_name = instrument.name
        else:
            instrument_name = str(instrument.__class__).replace("<class 'skorbit.instruments.", "").replace("'>", "")
    else:
        logging.error('Tried to get the instrument name from an object that was not a string or a subclass of PropagatorInstrument')
        raise TypeError
    return instrument_name


# why am I passing the name?
def _create_instrument(instrument: str | PropagatorInstrument, instrument_name: str,
                       tle_lines: List[str]) -> Optional[PropagatorInstrument]:
    """
    Handles the creation of default instruments. If passed an instance of
    the PropogatorInstrument class or it's subclasses, the instance is
    simply returned.

    instrument : str | PropagatorInstrument
        The instrument name as a string or a class instance of
        PropagatorInstrument or its subclasses. If a string only
        'aos-sky' and 'ticfire' are accepted.
    instrument_name: str
        The name of the instrument being created.
    tle_lines: List[str]
        The TLE as two strings in a list. Only used if creating a default
        instrument.

    Returns
    -------
    Optional[PropagatorInstrument]
        Returns the new instrument if one was created, the same instrument
        instance if one was passed in and None if an error occured.
    """
    inst = None
    # create the instrument
    if isinstance(instrument, str):
        if  instrument_name == 'aos-sky':
            inst = AOSSky(tle_lines)
            # print('Loaded AOS-Sky')
        elif instrument_name == 'ticfire':
            inst = TICFIRE(tle_lines)
            # print('Loaded TICFIRE')
        else:
            logging.error(f'Invlaid instrument: {instrument}. The only valid instruments for string arguments are "AOS-Sky" and "TICFIRE".')
            raise ValueError
    elif isinstance(instrument, PropagatorInstrument):
        inst = instrument
    else:
        logging.error('The instrument argument must be a string identifying an existing instrument '
                      'or an instance of a subclass of PropagatorInstrument')
        raise TypeError
    # if either error occurs then instrument will still be none
    return inst


def _convert_time(time: str | datetime) -> Optional[datetime]:
    """
    Converts the time passd in into a datetime object with a UTC timezone.

    Parameters
    ----------
    time : str | datetime
        Either a date and time as a string in the format
        'YYYY-MM-DD hh:mm:ss' or a datetime instance.

    Returns
    -------
    Optional[datetime]
        A datetime instance with the timezone set to UTC.
    """
    if isinstance(time, str):
        return read_time_string(time)
    elif isinstance(time, datetime):
        return time.replace(tzinfo=timezone.utc)
    else:
        logging.error('Tried to convert an object that was not a string or datetime into a datetime')
        raise TypeError


def _set_instrument_time(instrument: PropagatorInstrument, start: Optional[str | datetime], end: str | datetime, length: float, tle_time: datetime):
    """
    Configures the length of a PropagatorInstrument instance's simulation.

    Parameters
    ----------
    instrument : PropagatorInstrument
        The instrument to configure.
    start : Optional[str | datetime]
        The start time to configure the instrument with. If None then the
        epoch time of the TLE will be used. String times must have the
        format 'YYYY-MM-DD hh:mm:dd'.
    end : Optional[str | datetime]
        The end time to configure the instrument with. If None then the
        `length` keyword must be supplied. String times must have the
        format 'YYYY-MM-DD hh:mm:dd'.
    length : float
        The number of seconds from the start time to simulate the orbit for.
        This argument will be ignored if `end` is supplied. If this argument
        is None then `end` must be provided.
    tle_time : datetime
        The epoch time of the instrument TLE, as a datetime.
    """
    if start is not None:
        start_datetime = _convert_time(start)
    else:
        start_datetime = instrument.platform.platform_locator.time  # should be the time of the TLE

    if start_datetime is None:
        logging.error('Could not convert start argument into a datetime')
        raise TypeError

    if end is not None:
        end_datetime = _convert_time(end)
    else:
        end_datetime = start_datetime + timedelta(seconds=length)

    if end_datetime is None:
        logging.error('Could not convert end argument into a datetime')
        raise TypeError

    instrument.set_start_time((start_datetime - tle_time).total_seconds())
    instrument.set_end_time((end_datetime - tle_time).total_seconds())


def _make_filename_from_ds_time(instrument_name:str, ds: Dataset) -> str:
    """
    Generates and returns a name for an output file with the format:
    '{instrument_name}_{start_str}_to_{end_str}.nc'
    The start and end time are taken from the dataset. This is used
    for autonaming files when in multi-threaded mode.

    Parameters
    ----------
    instrument_name : str
        The name to use for the instrument.
    ds : Dataset
        A dataset with 'start_time' and 'end_time' attributes.

    Returns
    -------
    str
        The formatted name.
    """
    start_str = ds.start_time.replace(' ', '_').replace(':', '_')
    end_str = ds.end_time.replace(' ', '_').replace(':', '_')
    return f'{instrument_name}_{start_str}_to_{end_str}.nc'


def select_time_span(folder: Tuple[Path, str], start: datetime, end: datetime, search_pattern: str='*'):
    '''
    Searches through all files matching `search_pattern` and collects all data in the interval
    described by `start` and `end`. The collected data is returned as a new dataset, with the
    correct *start_time* and *end_time* attributes.

    Files from different instruments can only be combined if at least their *along* and
    *across* dimensions match. More generally, data from the instruments must be mergeable with
    `xarray.Dataset.merge()`.

    The function is unoptimized, and will load all files mathcing the pattern which can take a
    long time if there are many files.

    Parameters
    ----------
    folder : Tuple[Path, str]
        The folder to search for files in.
    start : datetime
        The start date and time of the interval to select.
    end : datetime
        The end date and time of the interval to select.
    search_pattern : str, optional
        The pattern to use when finding files to load. See `pathlib.Path.glob` for the rules of pattern
        matching.The default is '*', matching all files.

    Returns
    -------
    Dataset
        A dataset contianing all the data found for the interval from `start` to `end` in the files
        found in `folder`.
    '''
    if isinstance(folder, str):
        source_path = Path(folder)
    else:
        source_path = folder

    files = source_path.glob(search_pattern)

    datasets = []
    for file in files:
        ds = xr.load_dataset(file)
        selection = ds.sel(time=slice(start, end))
        if selection['time'].size > 0:
            datasets.append(selection)

    ds = xr.merge(datasets)
    ds.attrs['start_time'] = ds.time[0].dt.strftime('%Y-%m-%d %H:%M:%S').item()
    ds.attrs['end_time'] = ds.time[-1].dt.strftime('%Y-%m-%d %H:%M:%S').item()
    return ds


def compare_files(folder1: Tuple[Path, str], folder2: Tuple[Path, str], search_pattern: str ='*'):
    '''
    Compares the files in two folders for differences. The files in the second directory should
    have the same names as the files in the first directory, and should cover the same time period.

    The absolute difference is printed for variable in every pair of files that differs.

    Parameters
    ----------
    folder1 : Tuple[Path, str]
        The primary folder to load files from.
    folder2 : Tuple[Path, str]
        The secondary folder to load files from. Every file in folder1 that does not exist
        in this folder will be logged. Additional files in folder2 will not be noted.
    search_pattern : str, optional
        The pattern to use when loading files. See `pathlib.Path.glob` for the rules of pattern
        matching. The default is '*', matching all files.
    '''

    for file in folder1.glob(search_pattern):
        file2 = folder2.joinpath(file.name)
        if file2.exists():
            ds1 = xr.load_dataset(file)
            ds2 = xr.load_dataset(file2)
            if np.all(ds1['time'].values == ds2['time'].values):
                for key in ds1.keys():
                    difference = ds1[key].values.squeeze() - ds2[key].values.squeeze()
                    abs_diff = np.max(np.abs(difference)).item()
                    if abs_diff != 0:
                        print(f'{key} has an absolute difference of {abs_diff} for files: {file} {file2}')
            else:
                print(f'Times do not match for files {file} and {file2}')
        else:
            print(f'Pair for file: {file}, does not exist in: {folder2}')