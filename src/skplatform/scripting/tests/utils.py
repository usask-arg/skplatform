import numpy as np
from numpy.typing import NDArray


def print_all(obj):
    '''
    Quick and dirty debug command that prints all of an object's attributes
    '''
    for attr in dir(obj):
        try:
            print(attr, getattr(obj, attr))
        except:
            pass


def log10_histogram(data_array: NDArray, 
                    smallest_power: int = -18, largest_power: int = -3, num_bins=15,
                    header:bool = False, width: int = 3): #didn't work
    '''
    Calculates and prints a simple histogram, based on powers of 10:
    
    Example:
    
    print('Histogram:')
    position_histogram, histogram_header = log10_histogram(position_difference, header=True)
    velocity_histogram, _ = log10_histogram(velocity_difference)
    print(histogram_header, 'Log10')
    print(position_histogram, 'Position Relative Difference')
    print(velocity_histogram, 'Velocity Relative Difference')
    
    Output:

    Histogram:
    -18   -17   -16   -15   -14   -13   -12   -11   -10    -9    -8    -7    -6    -5    -4    -3    Log10
      |  0  |  0  |291  |288  | 80  |  7  |  2  |  0  |  0  |  0  |  0  |  0  |  0  |  0  |  0  | Position Relative Difference
      |  0  |  0  |266  |278  | 75  | 12  |  0  |  0  |  0  |  0  |  0  |  0  |  0  |  0  |  0  | Velocity Relative Difference
    '''
    # suppress annoying divide by zero and 0/0 errors
    with np.errstate(divide='ignore', invalid='ignore'):
        count, bins = np.histogram(np.log10(np.abs(data_array.flatten())), range=[smallest_power, largest_power], bins=num_bins)
        
        if header:
            header_txt = ''
            for k in range(len(bins)): 
                header_txt += f'{int(bins[k]):{width}}   '
        else: 
            header_txt = None

        data_txt = ' '*(width-1) + '|'
        for k in range(len(count)): 
            data_txt += f'{int(count[k]):{width}g}  |'
    
    return data_txt, header_txt