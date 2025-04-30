from typing import Optional
import click
from interface import simulate_orbit

# orbits break into files per orbit

@click.command()
@click.option('--instrument', default='aos-sky', type=str, help='The instrument to run. Must be either "AOS-Sky" or "TICFIRE". '
                                                                'The names are not case sensitive. The default is "aos-sky"')
@click.option('--output', default='', type=str, help='The path and filename to save the simulation output with. Can take on 3 '
                                                     'different values:\n\n'
                                                     'A directory: the file will be saved with an name using the format:\n'
                                                     '{instrument}_{starttime}_to_{endtime}.nc\n\n'
                                                     'An empty string: a file will be saved to the directory the script is running in, using automatic naming.\n\n'
                                                     'A full name and path: a file will be saved with give name and path.\n\n' 
                                                     'Any existing file will be overwritten.')
@click.option('--tle', default='', type=str, help='The file and path of a TLE to load. If not supplied the default TLE, '
                                                  'AOS_Sky_Descending_v2025-01-29_Epoch2019.tle, is used.')
@click.option('--start', default=None, type=str, help='The time to start the simulation from. This can be before the TLE epoch. '
                                                      'The format is strict, and must be "YYYY-MM-DD hh:mm:ss", quotes included. If not provided '
                                                      'the simulation starts with the epoch time in the TLE.')
@click.option('--end', default=None, type=str, help='The time to end the simulation at. This can be before the TLE epoch, '
                                                    'but only if --start is also provided. --start must be before --end. '
                                                    'The format is strict, and must be "YYYY-MM-DD hh:mm:ss", quotes included. '
                                                    'Overrides the --length argument.')
@click.option('--length', default=0, type=float, help='The number of seconds to run the simulation for. '
                                                      'Only used if --end is not set.') 
@click.option('--report_epoch', default=False, help='If set to 1, then the TLE epoch will be printed to the console and '
                                                    'the program will exit.')
@click.option('--threads', default=1, type=int, help='')

def run(instrument: str = '', output: str = '', tle: str = str, start: Optional[str] = None, 
        end: Optional[str] = None, 
        length: int = 0, report_epoch: bool = False, threads: int = 1):
    simulate_orbit(instrument=instrument, output=output, tle=tle, start=start, end=end,
                   length=length, report_epoch=report_epoch, threads=threads)


if __name__ == '__main__':
    run()