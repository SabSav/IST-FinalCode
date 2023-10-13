from SimulateRound import write_file
import os
import argparse
import numpy as np

if __name__ == '__main__':
    
  
    Ncells = 5
    r = 2.5
    compress_cell = 20
    dt = 1e-3
    sim_time = 100
    dir_path = os.path.dirname(os.path.realpath(__file__))
    sigma=0
    frame = 25
    A0 = 0.25
    h = r - np.sqrt(r**2 - A0*Ncells/np.pi)
   
    cases = [

      {'method': write_file, 'Ncells': Ncells, 'thick': h, 'r': r, 'dt': dt, 'sigma': sigma, 'apicaltension': 1.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':A0},
      {'method': write_file, 'Ncells': Ncells, 'thick': h, 'r': r, 'dt': dt, 'sigma': sigma, 'apicaltension': 2.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
        'dir_path': dir_path, 'frame_step': frame, 'A0':A0},
     {'method': write_file, 'Ncells': Ncells, 'thick': h, 'r': r, 'dt': dt, 'sigma': sigma, 'apicaltension': 3.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':A0},
     {'method': write_file, 'Ncells': Ncells, 'thick': h, 'r': r, 'dt': dt, 'sigma': sigma, 'apicaltension': 4.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':A0},
     {'method': write_file, 'Ncells': Ncells, 'thick': h, 'r': r, 'dt': dt, 'sigma': sigma, 'apicaltension': 5.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':A0},
     {'method': write_file, 'Ncells': Ncells, 'thick': h, 'r': r, 'dt': dt, 'sigma': sigma, 'apicaltension': 6.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':A0},
     {'method': write_file, 'Ncells': Ncells, 'thick': h, 'r': r, 'dt': dt, 'sigma': sigma, 'apicaltension': 7.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':A0},
     {'method': write_file, 'Ncells': Ncells, 'thick': h, 'r': r, 'dt': dt, 'sigma': sigma, 'apicaltension': 8.0, 'compress_cell': compress_cell, 'sim_time': sim_time, 
       'dir_path': dir_path, 'frame_step': frame, 'A0':A0},
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