from Advance import advanceLutolfExperiment
import numpy as np
import json
import argparse
from multiprocessing import Pool
from itertools import repeat

def initial_configuration(N_cells_cup, N_cells_flat, r_cup, L, h):
    cup_positive_x = []
    cup_positive_y = []
    cup_negative_x = []
    cup_negative_y = []
    flat = []
    for i in range(N_cells_cup + 1):
        if r_cup*np.cos(np.pi * i / (N_cells_cup)) >= 0:
            cup_positive_x.append(r_cup*np.cos(np.pi * i / (N_cells_cup)))
            cup_positive_y.append(r_cup*np.sin(np.pi * i / (N_cells_cup)))
        else:
            cup_negative_x.append(r_cup*np.cos(np.pi * i / (N_cells_cup)))
            cup_negative_y.append(r_cup*np.sin(np.pi * i / (N_cells_cup)))

    flat.append(L/N_cells_flat)
    for i in range(1, N_cells_flat-1):
        flat.append(flat[i-1] + L/N_cells_flat)
    
    basal = np.zeros((2*(len(cup_negative_x)+len(cup_positive_x)) + 2*len(flat), 2))
    #half cup positive
    cup_positive_x.sort()
    cup_positive_y.sort(reverse=True)
    basal[0:len(cup_positive_x), 0] = cup_positive_x
    basal[0:len(cup_positive_x), 1] = cup_positive_y
    #flat positive
    
    basal[len(cup_positive_x):len(cup_positive_x) + len(flat), 0] = r_cup
    basal[len(cup_positive_x):len(cup_positive_x) + len(flat), 1] = [-1*i for i in flat]

    #half cup positive traslated
    cup_positive_x.sort(reverse=True)
    cup_positive_y.sort()
    basal[len(cup_positive_x) + len(flat):2*len(cup_positive_x) + len(flat), 0] = cup_positive_x
    basal[len(cup_positive_x) + len(flat):2*len(cup_positive_x) + len(flat), 1] = [-1*(i+L) for i in cup_positive_y]

    #half cup negative traslated
    basal[2*len(cup_positive_x) + len(flat):2*len(cup_positive_x) + len(flat) + len(cup_negative_x), 0] = cup_negative_x
    basal[2*len(cup_positive_x) + len(flat):2*len(cup_positive_x) + len(flat) + len(cup_negative_x), 1] = [-1*(i+L) for i in cup_negative_y]
    
    flat.sort(reverse=True)
    #flat negative
    basal[2*len(cup_positive_x) + len(flat) + len(cup_negative_x): 2*len(cup_positive_x) + 2*len(flat) + len(cup_negative_x), 0] = -r_cup
    basal[2*len(cup_positive_x) + len(flat) + len(cup_negative_x): 2*len(cup_positive_x) + 2*len(flat) + len(cup_negative_x), 1] = [-1*i for i in flat]
    #half cup negative
    cup_negative_x.sort()
    cup_negative_y.sort()
    basal[2*len(cup_positive_x) + 2*len(flat) + len(cup_negative_x):2*len(cup_positive_x) + 2*len(flat) + 2*len(cup_negative_x), 0] = cup_negative_x
    basal[2*len(cup_positive_x) + 2*len(flat) + len(cup_negative_x):2*len(cup_positive_x) + 2*len(flat) + 2*len(cup_negative_x), 1] = cup_negative_y
            
    cup_positive_x = []
    cup_positive_y = []
    cup_negative_x = []
    cup_negative_y = []
    flat = []
    for i in range(N_cells_cup + 1):
        if (r_cup - h)*np.cos(np.pi * i / (N_cells_cup)) >= 0:
            cup_positive_x.append((r_cup - h)*np.cos(np.pi * i / (N_cells_cup)))
            cup_positive_y.append((r_cup - h)*np.sin(np.pi * i / (N_cells_cup)))
        else:
            cup_negative_x.append((r_cup - h)*np.cos(np.pi * i / (N_cells_cup)))
            cup_negative_y.append((r_cup - h)*np.sin(np.pi * i / (N_cells_cup)))

    flat.append(L/N_cells_flat)
    for i in range(1, N_cells_flat-1):
        flat.append(flat[i-1] + L/N_cells_flat)
    
    apical = np.zeros((2*(len(cup_negative_x)+len(cup_positive_x)) + 2*len(flat), 2))
    #half cup positive
    cup_positive_x.sort()
    cup_positive_y.sort(reverse=True)
    apical[0:len(cup_positive_x), 0] = cup_positive_x
    apical[0:len(cup_positive_x), 1] = cup_positive_y
    #flat positive
    apical[len(cup_positive_x):len(cup_positive_x) + len(flat), 0] = r_cup - h
    apical[len(cup_positive_x):len(cup_positive_x) + len(flat), 1] = [-1*i for i in flat]

    #half cup positive traslated
    cup_positive_x.sort(reverse=True)
    cup_positive_y.sort()
    apical[len(cup_positive_x) + len(flat):2*len(cup_positive_x) + len(flat), 0] = cup_positive_x
    apical[len(cup_positive_x) + len(flat):2*len(cup_positive_x) + len(flat), 1] = [-1*(i+L) for i in cup_positive_y]

    #half cup negative traslated
    apical[2*len(cup_positive_x) + len(flat):2*len(cup_positive_x) + len(flat) + len(cup_negative_x), 0] = cup_negative_x
    apical[2*len(cup_positive_x) + len(flat):2*len(cup_positive_x) + len(flat) + len(cup_negative_x), 1] = [-1*(i+L) for i in cup_negative_y]
    
    flat.sort(reverse=True)
    #flat negative
    apical[2*len(cup_positive_x) + len(flat) + len(cup_negative_x): 2*len(cup_positive_x) + 2*len(flat) + len(cup_negative_x), 0] = -(r_cup - h)
    apical[2*len(cup_positive_x) + len(flat) + len(cup_negative_x): 2*len(cup_positive_x) + 2*len(flat) + len(cup_negative_x), 1] = [-1*i for i in flat]
    #half cup negative
    cup_negative_x.sort()
    cup_negative_y.sort()
    apical[2*len(cup_positive_x) + 2*len(flat) + len(cup_negative_x):2*len(cup_positive_x) + 2*len(flat) + 2*len(cup_negative_x), 0] = cup_negative_x
    apical[2*len(cup_positive_x) + 2*len(flat) + len(cup_negative_x):2*len(cup_positive_x) + 2*len(flat) + 2*len(cup_negative_x), 1] = cup_negative_y

    return apical, basal
    
def simulate(Ncells, Ncellscup, Ncellsflat, L, rcup, h, apical_tension, compress_cell, dt, sigma, sim_time, frame, A0, seed):
    np.random.seed(seed)
    
    apical, basal= initial_configuration(N_cells_cup=Ncellscup, N_cells_flat=Ncellsflat, r_cup=rcup, L=L, h=h)
    config = advanceLutolfExperiment(N_cells=Ncells, apical=apical, basal=basal, N_cells_cup=Ncellscup, N_cells_flat=Ncellsflat, r_cup=rcup, L=L, h=h, apical_tension=apical_tension, compress_cell=compress_cell, dt=dt, sigma=sigma)    
    config.areas = config.evaluate_areas()
    for vertex in range(config.N_cells):
        config.A0[vertex] = A0
    bundle = config.export_dict()
    bundle['Number of cells on flat'] = config.N_cells_flat,
    bundle["Length"] = config.L
    bundle["Rcups"] = config.r_cup
    bundle["Thickness"] = config.h
    bundle['Target Area Flat'] = A0
    bundle['Target Area Curved'] = A0
    bundle['Apical Tension'] = config.apical_tension
    bundle['sigma'] = config.sigma
    energy = []
    config.energy()
    apical_time = np.zeros((int(sim_time//frame), config.N_cells, 2))
    basal_time = np.zeros((int(sim_time//frame), config.N_cells, 2))
    for t in np.arange(0, sim_time, config.dt):
        if t % frame == 0:
            index = int(t // frame) 
            apical_time[index, :, :] = config.apical
            basal_time[index, :, :] = config.basal
        energy.append(float(config.E))
        config.advance()
        print(f'Simulate Noise: Ncells={config.N_cells}, ApicalTension={config.apical_tension}, sigma={config.sigma}, seed={seed}, t = {t}, Energy = {config.E}')
    
    bundle['Energy'] = energy[:]
    bundle['Apical_x'] = apical_time[:, :, 0].tolist()
    bundle['Apical_y'] = apical_time[:, :, 1].tolist()
    bundle['Basal_x'] = basal_time[:, :, 0].tolist()
    bundle['Basal_y'] = basal_time[:, :, 1].tolist()
      
    return bundle

def merge(Ncells, Ncellscup, Ncellsflat, lflat, r, h, apicaltension, compress_cell, dt, sigma, sim_time, frame, A0, Nsamples):
    np.random.seed(0)
    seed = np.random.randint(2**32 - 1, size=Nsamples)
    with Pool(processes=16) as pool:
        simulation = pool.starmap(simulate, zip(repeat(Ncells), repeat(Ncellscup), repeat(Ncellsflat), repeat(lflat), repeat(r), repeat(h), repeat(apicaltension), 
                                                repeat(compress_cell), repeat(dt), repeat(sigma), repeat(sim_time), repeat(frame), repeat(A0), seed))
    return simulation

def write_file(Ncellscup, Ncellsflat, lflat, r, thick, apicaltension, compress_cell, dt, sigma, sim_time, frame, dir_path, A0, Nsamples):
        Ncells = 2*Ncellscup + 2*Ncellsflat
        bundle = merge(Ncells, Ncellscup, Ncellsflat, lflat, r, thick, apicaltension, compress_cell, dt, sigma, sim_time, frame, A0, Nsamples)
        for i in range(Nsamples):
            output = f'{dir_path}/Lutolf/Noise/DifferentDensity/{Ncells}/sigma/{sigma}/bundle-{i}.json'
            with open(output, 'w') as file:
                json.dump(bundle[i], file, indent=20)


if __name__ == '__main__':
    # Prepare arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dir_path', dest='dir_path', type=str, help="Directory for saving file")
    parser.add_argument(
        '-Ncellsflat', dest='Ncellsflat', type=int, default=3, help="Number of cells on flat side"
    )
    parser.add_argument(
        '-Ncellscup', dest='Ncellscup', type=int, default=3, help="Number of cells on curved side"
    )
    parser.add_argument(
        '-lfat', dest='L', type=float, default=0.0, help="Length flat side"
    )
    parser.add_argument(
        '-r', dest='r_cup', type=float, default=0.0, help="Radius cup (apical)"
    )
    parser.add_argument(
        '-thick', dest='h', type=float, default=0.0, help="Initial thickness"
    )
    parser.add_argument(
        '-dt', dest='dt', type=float, default=0.5,
        help="dt"
    )
    parser.add_argument(
        '-sigma', dest='sigma', type=float, default=0.5,
        help="noise"
    )
    parser.add_argument(
        '-apicaltension', dest='apicaltension', type=float, default=1.0,
        help="Apical Tension"
    )
    parser.add_argument(
        '-compress_cell', dest='compress_cell', type=int, default=10,
        help="Compressibility"
    )
    parser.add_argument(
        '-sim_time', dest='sim_time', type=int, default=10,
        help="Simulation Time"
    )
    parser.add_argument(
        '-Nsamples', dest='Nsamples', type=int, default=10,
        help="Number of samples"
    )
    parser.add_argument(
        '-frame', dest='frame', type=int, default=10,
        help="Frame"
    )
    parser.add_argument(
        '-A0', dest='A0', type=int, default=10,
        help="A0"
    )