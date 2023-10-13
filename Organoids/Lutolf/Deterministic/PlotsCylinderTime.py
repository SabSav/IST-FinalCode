import matplotlib.pyplot as plt
import json
import numpy as np
import os

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


if __name__ == '__main__':
    
    dir_path = os.path.dirname(os.path.realpath(__file__))
    Ncells=100
    apical_tension = np.arange(1.0, 8.5, 0.5)
    l=7.5
    rcup = 2.5
    Ncellscup = 25
    Ncellsflat = 25
    flat_lat = l/Ncellsflat
    h = 0.5/flat_lat

    apical, basal = initial_configuration(N_cells_cup=Ncellscup, N_cells_flat=Ncellsflat, r_cup=rcup, L=l, h=h)

    x = []
    y = []
    for i in range(Ncells):
        x.append(apical[i, 0])
        x.append(basal[i, 0])
        y.append(apical[i, 1])
        y.append(basal[i,1])
    for i in range(0, len(x),2):
        plt.plot(x[i:i+2], y[i:i+2], color="black", linewidth=0.5)
    plt.plot(rcup*np.cos(np.arange(0, np.pi, 1e-3)), rcup*np.sin(np.arange(0, np.pi, 1e-3)), color='darkred', label='basal')
    plt.plot(apical[:, 0], apical[:, 1], color='navy', label='apical')
    plt.plot([apical[0, 0], apical[-1, 0]], [apical[0, 1], apical[-1, 1]], color="navy")
    plt.vlines(rcup, -l, 0, color='darkred')
    plt.vlines(-rcup, -l, 0, color='darkred')
    plt.plot(rcup*np.cos(np.arange(0, np.pi, 1e-3)), -(rcup*np.sin(np.arange(0, np.pi, 1e-3)) + l), color='darkred')
    plt.xlim(-10,10)
    plt.ylim(-10.5, 3)
    plt.axis('off')
    plt.legend()
    plt.title(f'Initial configuration')  
    
    for at in apical_tension:
        plt.savefig(f'{dir_path}/{Ncells}/ApicalTension/{at}/InitialConfig.jpg')

    plt.close()    

    for at in apical_tension:
        file = open(f'{dir_path}/{Ncells}/ApicalTension/{at}/bundle.json')
        data = json.load(file)
        Ncells = data['Number of cells']
        apical_y = data['Apical_y'][-1]
        apical_x = data['Apical_x'][-1]
        basal_x = data['Basal_x'][-1]
        basal_y = data['Basal_y'][-1]
        engy = data['Energy']
        
        x = []
        y = []
        for i in range(Ncells):
            x.append(apical_x[i])
            x.append(basal_x[i])
            y.append(apical_y[i])
            y.append(basal_y[i])
        for i in range(0, len(x),2):
            plt.plot(x[i:i+2], y[i:i+2], color="black", linewidth=0.5)
        plt.plot(rcup*np.cos(np.arange(0, np.pi, 1e-3)), rcup*np.sin(np.arange(0, np.pi, 1e-3)), color='darkred', label='basal')
        plt.plot(apical_x[:], apical_y[:], color='navy', label='apical')
        plt.plot([apical_x[0], apical_x[-1]], [apical_y[0], apical_y[-1]], color="navy")
        plt.vlines(rcup, -l, 0, color='darkred')
        plt.vlines(-rcup, -l, 0, color='darkred')
        plt.plot(rcup*np.cos(np.arange(0, np.pi, 1e-3)), -(rcup*np.sin(np.arange(0, np.pi, 1e-3)) + l), color='darkred')
        plt.xlim(-10,10)
        plt.ylim(-10.5, 3)
        plt.axis('off')
        plt.legend()
        plt.title(f'Final Configuration') 
        plt.savefig(f'{dir_path}/{Ncells}/ApicalTension/{at}/FinalConfig.jpg')
        plt.close()

        plt.plot(engy)
        plt.savefig(f'{dir_path}/{Ncells}/ApicalTension/{at}/Energy.jpg')
