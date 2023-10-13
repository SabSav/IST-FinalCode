import matplotlib.pyplot as plt
import json
import numpy as np
import math
import os
import seaborn as sns
import pandas as pd


def fhc(x, Af, Ac, L, R, N, at):
    return ((np.pi**2)*Af/(2*L*Ac**2))*x**4 - 2*np.pi*R/Ac*(0.5 + np.pi*Af/(Ac*L))*x**3 + 2*np.pi/Ac*(R + (Af/L)*(0.5*N + np.pi*R**2/Ac))*x**2 - 2*np.pi*(N*Af*R/(Ac*L) + at)*x
if __name__ == '__main__':

    apical_tension = np.arange(1.0, 4.0, 0.5)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    Ncells=20
    for at in apical_tension:
        file = open(f'{dir_path}/DifferentDensity/{Ncells}/ApicalTension/{at}/bundle.json')
        data = json.load(file)
        Ncellsflat = data['Number of cells on flat']
        rcup = data['Rcups']
        apical_y = data['Apical_y'][-1]
        apical_x = data['Apical_x'][-1]
        basal_x = data['Basal_x'][-1]
        basal_y = data['Basal_y'][-1]
        l = data['Length']
        h = data['Thickness']
        hcurvedSimulation = []
        lapical = []
        y_flat_plus = []
        y_flat_neg = []
        x_curved_plus = []
        x_curved_neg = []
        count_curved = 0
        for vertex in range(Ncells):
            if basal_x[vertex] == rcup and -l < basal_y[vertex] < 0 :
                y_flat_plus.append(basal_y[vertex])
            elif basal_x[vertex] == -rcup and -l < basal_y[vertex] < 0:
                y_flat_neg.append(basal_y[vertex])
            else:
                if basal_y[vertex] >= 0:
                    x_curved_plus.append(basal_x[vertex])
                elif basal_y[vertex] <= -l:
                    x_curved_neg.append(basal_x[vertex])
        y_flat_cross_plus = [max(y_flat_plus), min(y_flat_plus)]
        y_flat_cross_neg = [max(y_flat_neg), min(y_flat_neg)]
        x_curved_cross_neg = [max(x_curved_neg), min(x_curved_neg)]
        x_curved_cross_plus = [max(x_curved_plus), min(x_curved_plus)]
        for vertex in range(Ncells):
            if (basal_x[vertex] == rcup and basal_y[vertex] not in y_flat_cross_plus) or (basal_x[vertex] == -rcup and basal_y[vertex] not in y_flat_cross_neg):
                lapical.append(np.sqrt((apical_x[vertex] - apical_x[(vertex +1)%Ncells])**2 + (apical_y[vertex] - apical_y[(vertex +1)%Ncells])**2))
            elif (basal_y[vertex] >= 0 and basal_x[vertex] not in x_curved_cross_plus) or (basal_y[vertex] <= -l and basal_x[vertex] not in x_curved_cross_neg):
                    hcurvedSimulation.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
                    count_curved += 1
        lcurved = []
        for i in range(len(hcurvedSimulation)):
            lcurved.append(np.pi*(rcup-0.5*hcurvedSimulation[i])/(0.5*count_curved))
        lcurved = np.mean(lcurved)*10
        lapical = np.max(lapical)*10
        lapical_tot = [lcurved, lapical]
        regions = ['Tips', 'Sides']
        data = {'Internuclear Distance':lapical_tot, 'Region': regions}
        df = pd.DataFrame(data)
        sns.catplot(x='Region', y="Internuclear Distance", data=df, color = "black", s=50)
        plt.ylim(10, 50)
        plt.savefig(f'{dir_path}/DifferentDensity/{Ncells}/ApicalTension/{at}/InternuclearDistance.jpg')



    
    