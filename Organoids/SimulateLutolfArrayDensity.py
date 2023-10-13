from SimulateLutolfDifferentDensity import write_file
import os
import argparse
import numpy as np

if __name__ == '__main__':
    
  
    Ncellscup = 9
    Ncellsflat = 9
    A0 = 0.5
    lflat = 7.5
    flat_lat = lflat/Ncellsflat
    h = A0/flat_lat
    r = 2.5
    compress_cell = 20
    dt = 1e-3
    sim_time = 500
    dir_path = os.path.dirname(os.path.realpath(__file__))
    sigma=0
    frame = 25
    
    
    cases = [
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 1.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 1.5, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 2.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 2.5, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 3.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 3.5, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 4.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 4.5, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 5.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 5.5, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 6.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 6.5, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 7.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 7.5, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5},
     {'method': write_file, 'Ncellsflat': Ncellsflat, 'Ncellscup': Ncellscup, 'thick': h, 'r': r, 'lflat': lflat,
      'dt': dt, 'sigma': sigma, 'apicaltension': 8.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':0.5}
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