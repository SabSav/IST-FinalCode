from Advance import advanceRound
import numpy as np
import json
import argparse

def initial_configuration(N_cells, r, h):
        
    basal = np.zeros((N_cells, 2))
    apical = np.zeros((N_cells, 2))

    for i in range(N_cells):
        basal[i, 0] = r*np.cos(-2*np.pi * i / N_cells)
        basal[i, 1] = r*np.sin(-2*np.pi * i / N_cells)
        apical[i, 0] = (r-h)*np.cos(-2*np.pi * i / N_cells)
        apical[i, 1] = (r-h)*np.sin(-2*np.pi * i / N_cells)

    return apical, basal



def simulate(N_cells, r, h, A0, apical_tension, compress_cell, dt, sigma, sim_time, frame_step):
    
    apical, basal =initial_configuration(N_cells, r, h)
    config = advanceRound(N_cells=N_cells, apical=apical, basal=basal, r=r, h=h, apical_tension=apical_tension, compress_cell=compress_cell, dt=dt, sigma=sigma)
    for vertex in range(config.N_cells):
         config.A0[vertex] = A0
    config.areas = config.evaluate_areas()
    lapical = np.zeros(config.N_cells)
    for vertex in range(config.N_cells):
        lapical[vertex] = np.sqrt((config.apical[vertex, 0] - config.apical[(vertex + 1)%config.N_cells, 0])**2 + (config.apical[vertex, 1] - config.apical[(vertex + 1)%config.N_cells, 1])**2)
    
    bundle = config.export_dict()
    bundle["Radius"] = config.r
    bundle["Thickness"] = config.h
    bundle['Target Area'] = config.A0
    bundle['Initial areas'] = config.areas.tolist()
    bundle['Apical Lenght init'] = lapical.tolist()
    bundle['Apical Tension'] = config.apical_tension
    bundle['sigma'] = config.sigma
    energy = []
    apical_time = np.zeros((int(sim_time//frame_step), config.N_cells, 2))
    basal_time = np.zeros((int(sim_time//frame_step), config.N_cells, 2))
    
    for t in np.arange(0, sim_time, config.dt):
        if t % frame_step == 0:
            index = int(t // frame_step) 
            apical_time[index, :, :] = config.apical
            basal_time[index, :, :] = config.basal
        config.advance()
        energy.append(float(config.E))
        print(f'Simulate Round Organoid, r={config.r}, Ncells={config.N_cells}, A0 = {config.A0[0]}, K={config.compress_cell}, at={config.apical_tension}, t = {t} , Energy = {config.E}')
    lapical = np.zeros(config.N_cells)
    for vertex in range(config.N_cells):
        lapical[vertex] = np.sqrt((config.apical[vertex, 0] - config.apical[(vertex + 1)%config.N_cells, 0])**2 + (config.apical[vertex, 1] - config.apical[(vertex + 1)%config.N_cells, 1])**2)
    config.areas = config.evaluate_areas()
    bundle['Final areas'] = config.areas.tolist()
    bundle['Apical Lenght final'] = lapical.tolist()
    bundle['Energy'] = energy[:]
    bundle['Apical_x'] = apical_time[:, :, 0].tolist()
    bundle['Apical_y'] = apical_time[:, :, 1].tolist()
    bundle['Basal_x'] = basal_time[:, :, 0].tolist()
    bundle['Basal_y'] = basal_time[:, :, 1].tolist()
      
    return bundle
    

def write_file(Ncells, r, thick, A0, apicaltension, compress_cell, dt, sigma, sim_time, dir_path, frame_step):
        bundle = simulate(Ncells, r, thick, A0, apicaltension, compress_cell, dt, sigma, sim_time, frame_step)
        output = f'{dir_path}/Round/DifferentRadius/{r}/DifferentNcells/{Ncells}/ApicalTension/{apicaltension}/bundle_A0={A0}.json'
        with open(output, 'w') as file:
            json.dump(bundle, file, indent=20)


if __name__ == '__main__':
    # Prepare arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dir_path', dest='dir_path', type=str, help="Directory for saving file")
    parser.add_argument(
        '-Ncells', dest='Ncells', type=int, default=3, help="Number of cells"
    )
    parser.add_argument(
        '-r', dest='r', type=float, default=0.0, help="Radius"
    )
    parser.add_argument(
        '-thick', dest='h', type=float, default=0.0, help="Initial thickness"
    )
    parser.add_argument(
        '-A0', dest='A0', type=float, default=0.0, help="Preferred Area"
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
        '-frame_step', dest='frame_step', type=int, default=10,
        help="Frame Step"
    )