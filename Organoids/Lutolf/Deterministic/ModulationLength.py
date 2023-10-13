import matplotlib.pyplot as plt
import json
import numpy as np
import math
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import os
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import seaborn as sns

def fhc(x, Af, Ac, L, R, N, at):
    return ((np.pi**2)*Af/(2*L*Ac**2))*x**4 - 2*np.pi/Ac*(0.5 + R*np.pi*Af/(Ac*L))*x**3 + 2*np.pi/Ac*(R + (Af/L)*(0.5*N + np.pi*R**2/Ac))*x**2 - 2*np.pi*(N*Af*R/(Ac*L) + at)*x
if __name__ == '__main__':

    apical_tension = np.arange(1.0, 3.5, 0.5)
    lenghts = [3.75, 7.5, 10, 15, 20]
    dir_path = os.path.dirname(os.path.realpath(__file__))
    hratioMean = np.zeros((len(apical_tension), len(lenghts)))
    hratioStd = np.zeros((len(apical_tension), len(lenghts)))
    mod = np.zeros((len(apical_tension), len(lenghts)))
    hc = np.zeros((len(apical_tension), len(lenghts)))
    hf = np.zeros((len(apical_tension), len(lenghts)))
    rcup = 2.5
    A_flat = 0.5
    A_curved = 0.5
    for at in range(len(apical_tension)):
        for l in range(len(lenghts)):
            file = open(f'{dir_path}/DifferentShape/DifferentLength/{lenghts[l]}/ApicalTension/{apical_tension[at]}/bundle.json')
            data = json.load(file)
            Ncells = data["Number of cells"]
            Ncellsflat = data['Number of cells on flat']
            apical_y = data['Apical_y'][-1]
            apical_x = data['Apical_x'][-1]
            basal_x = data['Basal_x'][-1]
            basal_y = data['Basal_y'][-1]
            h = data['Thickness']
            hcross = []
            hcurvedSimulation = []
            hflatSimulation = []
            y_flat_plus = []
            y_flat_neg = []
            x_curved_plus = []
            x_curved_neg = []        
            for vertex in range(Ncells):
                if basal_x[vertex] == rcup and -lenghts[l] < basal_y[vertex] < 0 :
                    y_flat_plus.append(basal_y[vertex])
                elif basal_x[vertex] == -rcup and -lenghts[l] < basal_y[vertex] < 0:
                    y_flat_neg.append(basal_y[vertex])
                else:
                    if basal_y[vertex] >= 0:
                        x_curved_plus.append(basal_x[vertex])
                    elif basal_y[vertex] <= -lenghts[l]:
                        x_curved_neg.append(basal_x[vertex])
            y_flat_cross_plus = [max(y_flat_plus), min(y_flat_plus)]
            y_flat_cross_neg = [max(y_flat_neg), min(y_flat_neg)]
            x_curved_cross_neg = [max(x_curved_neg), min(x_curved_neg)]
            x_curved_cross_plus = [max(x_curved_plus), min(x_curved_plus)]
            for vertex in range(Ncells):
                if len(y_flat_neg) > 2 and len(y_flat_plus) > 2:
                    if basal_y[vertex] in y_flat_cross_plus or basal_y[vertex] in y_flat_cross_neg:
                        hcross.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
                    elif basal_x[vertex] in x_curved_cross_plus or basal_x[vertex] in x_curved_cross_neg:
                        hcross.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
                    elif  (basal_x[vertex] == rcup and basal_y[vertex] not in y_flat_cross_plus) or (basal_x[vertex] == -rcup and basal_y[vertex] not in y_flat_cross_neg):
                        hflatSimulation.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
                    else:
                        hcurvedSimulation.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
                else:
                    if basal_x[vertex] == rcup or basal_x[vertex] == -rcup:
                        hflatSimulation.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
                    elif basal_x[vertex] in x_curved_cross_neg or basal_x[vertex] in x_curved_cross_plus:
                        basal_x[vertex] = basal_x[vertex]
                    else:
                        hcurvedSimulation.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
            
            hratioMean[at, l] = np.mean(hflatSimulation)*np.mean([i**(-1) for i in hcurvedSimulation])
            hratioStd[at, l] = np.mean([i**2 for i in hflatSimulation])*np.mean([(i**(-1))**2 for i in hcurvedSimulation]) - (np.mean(hflatSimulation)**2)*(np.mean([i**(-1) for i in hcurvedSimulation]))**2
            mod[at, l] = -hratioMean[at, l] + 1
            res = minimize(fhc, x0=h, args=(A_flat, A_curved, lenghts[l], rcup, Ncells, apical_tension[at]))
            hc[at, l] = res.x
            hf[at, l] = 0.5*Ncells*A_flat/lenghts[l] - (A_flat*np.pi/(lenghts[l]*A_curved))*(rcup*hc[at, l] - 0.5*hc[at, l]**2)
    
    sns.set_style("whitegrid")
    plt.figure(figsize=(8,7))
    color = []
    legends = []
    np.random.seed(0)
    for i in range(len(apical_tension)):
        color.append('#%06X' % np.random.randint(0, 0xFFFFFF))
        legends.append(mpatches.Patch(color=color[i], label=f'Γa = {apical_tension[i]}'))
    for at in range(len(apical_tension)):
        X_r_hratio = interp1d(lenghts, hf[at, :]/hc[at, :], kind='cubic')
        X_r = np.linspace(min(lenghts), max(lenghts), 500)
        hratio = X_r_hratio(X_r)
        plt.plot(X_r, -hratio+1, color=color[at])
        plt.errorbar(lenghts, mod[at, :], yerr=hratioStd[at, :], fmt='o', color=color[at], label = f'Simulation')
    line = Line2D([0], [0], label='Theory', color='k')
    point = Line2D([0], [0], color='k', marker='o', markersize=5,ls="", label='Simulation')
    line_exp = Line2D([0], [0], label='Experiment', color='gray', linestyle='dashed')
    legends.append(line)
    legends.append(point)
    legends.append(line_exp)
    plt.hlines(y=0.6, xmin=min(lenghts), xmax=16, linestyles='dashed', colors='gray')
    plt.hlines(y=0.7, xmin=min(lenghts), xmax=16, linestyles='dashed', colors='gray')
    plt.legend(handles=legends, loc='upper right', framealpha=0.1)
    plt.xlabel('Length')
    plt.ylabel(r'$\Omega$', fontsize=15)
    plt.savefig(f'{dir_path}/DifferentShape/DifferentLength/ModLength.jpg')