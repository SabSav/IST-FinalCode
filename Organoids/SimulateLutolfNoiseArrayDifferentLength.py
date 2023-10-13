from SimulateLutolfNoiseDifferentLength import write_file
import os
import argparse
import numpy as np

if __name__ == '__main__':
    
    r = 2.5
    lflat = 10
    lflat0 = 7.5
    density = 20/(2*np.pi*r + 2*lflat0) #has to be constant
    Ncellsflat = round(density*lflat)
    Ncellscup = round(density*np.pi*r)
    A0 = 0.5
    h = A0/(lflat/Ncellsflat)
    compress_cell = 20
    sim_time = 500
    dir_path = os.path.dirname(os.path.realpath(__file__))
    frame = 250
    dt = 1e-3
    apicaltension=3.5

    Nsamples = 100
    
    cases = [
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': 0.01, 'apicaltension': apicaltension, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'Nsamples': Nsamples, 'frame': frame, 'A0':A0},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': 0.05, 'apicaltension': apicaltension, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'Nsamples': Nsamples, 'frame': frame,'A0':A0},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': 1e-5, 'sigma': 0.1, 'apicaltension': apicaltension, 'compress_cell': compress_cell, 'sim_time': 20, 
       'dir_path': dir_path, 'Nsamples': Nsamples, 'frame': 10,'A0':A0},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': 1e-4, 'sigma': 0.2, 'apicaltension': apicaltension, 'compress_cell': compress_cell, 'sim_time': 100, 
       'dir_path': dir_path, 'Nsamples': Nsamples, 'frame': 50, 'A0':A0},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': 1e-4, 'sigma': 0.5, 'apicaltension': apicaltension, 'compress_cell': compress_cell, 'sim_time': 100, 
       'dir_path': dir_path, 'Nsamples': Nsamples, 'frame': 50, 'A0':A0}
    ]

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '-i', dest='cases', type=str, default=':',
        help="Case indices (default: ':')"
    )
    args = parser.parse_args()

    for case in eval(f'np.asarray(cases)[{args.cases}]'):
        kwargs = case.copy()
        method = kwargs.pop('method')
        method(**kwargs)
        print()
