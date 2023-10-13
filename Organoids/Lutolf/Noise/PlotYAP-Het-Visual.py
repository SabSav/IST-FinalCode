import json
import numpy as np
import os
import math
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon


if __name__ == '__main__':


    dir_path = os.path.dirname(os.path.realpath(__file__))
    l = 7.5
    Ncells = 20
    A0 = 0.5
    sigma = 0.05
    rcup = 2.5
    Nsample = 100
    time = 500
    dt = 1e-3
    density = Ncells/(2*np.pi*rcup + 2*l)
    yap = np.zeros((Nsample, Ncells))
    het = np.zeros((Nsample, Ncells))
    engy = np.zeros((Nsample, int(time/dt)))
    apical_x_mean = np.zeros((Nsample, Ncells))
    apical_y_mean = np.zeros((Nsample, Ncells))
    basal_x_mean = np.zeros((Nsample, Ncells))
    basal_y_mean = np.zeros((Nsample, Ncells))
    for j in range(Nsample):
        file = open(f'{dir_path}/DifferentDensity/{Ncells}/sigma/{sigma}/bundle-{j}.json')
        data = json.load(file)
        apical_y_mean[j, :] = data['Apical_y'][-1]
        apical_x_mean[j, :] = data['Apical_x'][-1]
        basal_x_mean[j, :]= data['Basal_x'][-1]
        basal_y_mean[j, :] = data['Basal_y'][-1]
        engy[j] = data['Energy']

        for vertex in range(Ncells):
            la = math.sqrt((apical_x_mean[j, vertex] - apical_x_mean[j, (vertex + 1)%Ncells])**2 + (apical_y_mean[j, vertex] - apical_y_mean[j, (vertex + 1)%Ncells])**2)
            yap[j, vertex] = la
        for vertex in range(Ncells):
            het[j, vertex] = 0.5*((yap[j, vertex] - yap[j, (vertex+1)%Ncells])**2 + (yap[j, vertex] - yap[j, (vertex-1)%Ncells])**2)
    basal_x_mean = np.mean(basal_x_mean, axis=0)
    basal_y_mean = np.mean(basal_y_mean, axis=0)
    apical_x_mean = np.mean(apical_x_mean, axis=0)
    apical_y_mean = np.mean(apical_y_mean, axis=0)
    yap = np.mean(yap, axis=0)
    yap = [1/(max(yap)- min(yap))*(i-min(yap)) for i in yap]
    het = np.mean(het, axis=0)

    plt.figure()
    engy = np.mean(engy, axis=0)
    plt.plot(engy)
    plt.savefig(f'{dir_path}/DifferentDensity/{Ncells}/sigma/{sigma}/Energy.jpg')
    plt.close()
    
    centroid = np.zeros((Ncells, 2))
    for vertex in range(Ncells):
        x = [basal_x_mean[vertex], basal_x_mean[(vertex + 1)%Ncells], apical_x_mean[vertex], apical_x_mean[(vertex + 1)%Ncells]]
        y = [basal_y_mean[vertex], basal_y_mean[(vertex + 1)%Ncells], apical_y_mean[vertex], apical_y_mean[(vertex + 1)%Ncells]]
        centroid[vertex, 0] = sum(x) / 4
        centroid[vertex, 1] = sum(y) / 4
    fig, ax = plt.subplots(1)
    a = ax.scatter(centroid[:, 0], centroid[:, 1], c=yap, cmap='magma')
    cba = fig.colorbar(a)
    x = []
    y = []
    for i in range(Ncells):
        x.append(apical_x_mean[i])
        x.append(basal_x_mean[i])
        y.append(apical_y_mean[i])
        y.append(basal_y_mean[i])
    for i in range(0, len(x),2):
        ax.plot(x[i:i+2], y[i:i+2], color="black", linewidth=0.5)
    a1 = ax.plot(basal_x_mean[0:Ncells], basal_y_mean[0:Ncells], color="black", linewidth=0.5)
    a1 = ax.plot([basal_x_mean[0], basal_x_mean[-1]], [basal_y_mean[0], basal_y_mean[-1]], color="black", linewidth=0.5)
    a1 = ax.plot(apical_x_mean[0:Ncells], apical_y_mean[0:Ncells], color="black", linewidth=0.5)
    a1 = ax.plot([apical_x_mean[0], apical_x_mean[-1]], [apical_y_mean[0], apical_y_mean[-1]], color="black", linewidth=0.5)
    a1 = ax.set_xlim(-10,10)
    a1 = ax.set_ylim(-(l+rcup+0.5),rcup+0.5)
    a = ax.axis('off')
    a1 = ax.set_title(f'YAP level, Density = {round(Ncells/(2*np.pi*rcup + 2*l), 2)}, σ = {sigma} (over {Nsample} organoids)')
    fig.savefig(f'{dir_path}/DifferentDensity/{Ncells}/sigma/{sigma}/YAPnuclei.jpg')
    


    fig, ax = plt.subplots(1)
    patches = []
    for i in range(0, len(x),2):
        ax.plot(x[i:i+2], y[i:i+2], color="black", linewidth=0.5)
    ax.plot(basal_x_mean[0:Ncells], basal_y_mean[0:Ncells], color="black", linewidth=0.5)
    ax.plot([basal_x_mean[0], basal_x_mean[-1]], [basal_y_mean[0], basal_y_mean[-1]], color="black", linewidth=0.5)
    ax.plot(apical_x_mean[0:Ncells], apical_y_mean[0:Ncells], color="black", linewidth=0.5)
    ax.plot([apical_x_mean[0], apical_x_mean[-1]], [apical_y_mean[0], apical_y_mean[-1]], color="black", linewidth=0.5)
    for vertex in range(Ncells):
        verts = [(apical_x_mean[vertex], apical_y_mean[vertex]), (basal_x_mean[vertex], basal_y_mean[vertex]), (basal_x_mean[(vertex+1) % Ncells], basal_y_mean[(vertex + 1) % Ncells]), (apical_x_mean[(vertex+1) % Ncells], apical_y_mean[(vertex + 1) % Ncells])]
        polygon = Polygon(verts, closed=True)
        patches.append(polygon)
    collection = PatchCollection(patches, cmap = 'Reds', match_original=True)
    collection.set_array(het)
    a = ax.add_collection(collection)
    fig.colorbar(a)
    a = ax.set_xlim(-10, 10)
    a = ax.set_ylim(-(l+rcup + 0.5), rcup + 0.5)
    a = ax.axis('off')
    a = ax.set_title(f'YAP Heterogeneity (over {Nsample} organoids),  σ={sigma} ')
    fig.savefig(f'{dir_path}/DifferentDensity/{Ncells}/sigma/{sigma}/HetVisual.jpg')

    