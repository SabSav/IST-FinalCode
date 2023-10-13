import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import os
import json
import numpy as np
import math 
import seaborn as sns

if __name__ == '__main__':  
    
    def sigmoid(x, a,b,c,d):
        return c + d*(1/ (1 + np.exp(-x/b + a)))
    def derivative(x, a,b,d):
        return (d/b)*(np.exp(-x/b + a)/((1 + np.exp(-x/b + a)))**2)
    

    dir_path = os.path.dirname(os.path.realpath(__file__))
    l = 7.5
    sigma = 0.05
    Nsamples = 100
    A0 = 0.5
    dt = 1e-3
    z = []
    rcup = 2.5
    file = open(f'{dir_path}/A0cases/{A0}/DifferentShape/DifferentLength/{l}/sigma/{sigma}/bundle-{0}.json')
    data = json.load(file)
    Ncells = data['Number of cells']
    Nf = data['Number of cells on flat'] + 1
    Nc = int(0.5*(0.5*Ncells - Nf))
    theta = np.linspace(0, 0.5*np.pi, Nc)
    for i in range(len(theta)):
        z.append(np.cos(theta[i])*rcup)
    vertex = len(z) - 1
    for i in range(vertex + 1, vertex + 1 + Nf):
        z.append(z[i - 1] - l/Nf)
    theta = np.flip(theta)
    for i in range(1, len(theta)):
        z.append(-l - np.cos(theta[i])*rcup)    
    z = np.array(z)
    z = np.array([(i - max(z))/(min(z) - max(z)) for i in z])
    util = []
    
    for i in range(len(z)):
        if 0 <= z[i] <= 0.5:
            util.append(i)
    sns.set_style("whitegrid")

    lapical = np.zeros((Nsamples, Ncells))
    lab_diff = np.zeros((Nsamples, Ncells))
    lab_ratio = np.zeros((Nsamples, Ncells))
    for j in range(Nsamples):
        file = open(f'{dir_path}/A0cases/{A0}/DifferentShape/DifferentLength/{l}/sigma/{sigma}/bundle-{j}.json')
        data = json.load(file)  
        apical_y = data['Apical_y'][-1]
        apical_x = data['Apical_x'][-1]
        basal_x = data['Basal_x'][-1]
        basal_y = data['Basal_y'][-1]
        
        for vertex in range(Ncells):
            la = math.sqrt((apical_x[vertex] - apical_x[(vertex+1)%Ncells])**2 + (apical_y[vertex] - apical_y[(vertex+1)%Ncells])**2)
            lb = math.sqrt((basal_x[vertex] - basal_x[(vertex+1)%Ncells])**2 + (basal_y[vertex] - basal_y[(vertex+1)%Ncells])**2)
            lapical[j, vertex] = la
            lab_diff[j, vertex] = lb - la
            lab_ratio[j, vertex] = la/lb
    
    lapicalmean = np.mean(lapical, axis=0)
    std_lapical = np.std(lapical, axis=0)


    fig, ax = plt.subplots(1, 3, figsize=(20, 5))
    pars, cov = curve_fit(f=sigmoid, xdata=z[util], ydata=lapicalmean[0:len(z[util])], bounds=(-np.inf, np.inf), maxfev=10**8)
    Y_lapical = sigmoid(np.arange(min(z[util]), max(z[util]), 10e-4), *pars)
    X_theta = np.linspace(z[util].min(), z[util].max(), len(Y_lapical))
    X_Y_cubic_lapical = interp1d(z[util], std_lapical[0:len(z[util])], kind='linear')
    Y_std_lapical = X_Y_cubic_lapical(X_theta)
    ax[0].plot(np.arange(min(z[util]), max(z[util]), 10e-4), Y_lapical, color='navy')
    ax[0].errorbar(z[util], lapicalmean[0:len(z[util])], yerr=std_lapical[0:len(z[util])], fmt='o', color='navy')
    ax[0].fill_between(X_theta, Y_lapical - Y_std_lapical, Y_lapical + Y_std_lapical, facecolor='navy', alpha=0.3)
    ax[0].set_ylabel(r'$\overline{l_{apical}}$', fontsize=15)
    ax[0].set_xlabel('Position (z/L)', fontsize=15)
    
    ax[1].plot(Y_lapical, Y_std_lapical, color='navy')
    ax[1].scatter(lapicalmean, std_lapical, color='navy')
    ax[1].set_ylabel(r'$\sigma_{l_{apical}}$', fontsize=15)
    ax[1].set_xlabel(r'$\overline{l_{apical}}$', fontsize=15)
    
    der = derivative(np.linspace(z[util].min(), z[util].max(), len(Y_lapical)), pars[0], pars[1], pars[3])
    Y_stderror_lapical = Y_std_lapical/der
    der = derivative(z[util], pars[0], pars[1], pars[3])
    std_error_lapical = std_lapical[0:len(z[util])]/der
    ax[2].plot(X_theta, Y_stderror_lapical, color='navy')
    ax[2].scatter(z[util], std_error_lapical, color='navy')
    ax[2].set_ylabel(r'$\sigma_{z}(z)$', fontsize=15)
    ax[2].set_xlabel('Position (z/L)', fontsize=15)
    ax[2].hlines(y=0.05, xmin=0, xmax=0.5, linestyle='dashed', color='gray')

    fig.suptitle(f'Length flat side = {l}, σ = {sigma} (over 100 samples)')
    fig.savefig(f'{dir_path}/A0cases/{A0}/DifferentShape/DifferentLength/{l}/sigma/{sigma}/PI-La.jpg', bbox_inches='tight')

    lab_diff_mean = np.mean(lab_diff, axis=0)
    std_lab_diff = np.std(lab_diff, axis=0)
    fig, ax = plt.subplots(1, 3, figsize=(20, 5))
    pars, cov = curve_fit(f=sigmoid, xdata=z[util], ydata=lab_diff_mean[0:len(z[util])], bounds=(-np.inf, np.inf), maxfev=10**8)
    Y_lab_diff = sigmoid(np.arange(min(z[util]), max(z[util]), 10e-4), *pars)
    X_theta = np.linspace(z[util].min(), z[util].max(), len(Y_lab_diff))
    X_Y_cubic_lab_diff = interp1d(z[util], std_lab_diff[0:len(z[util])], kind='linear')
    Y_std_lab_diff = X_Y_cubic_lab_diff(X_theta)
    ax[0].plot(np.arange(min(z[util]), max(z[util]), 10e-4), Y_lab_diff, color='darkred')
    ax[0].errorbar(z[util], lab_diff_mean[0:len(z[util])], yerr=std_lab_diff[0:len(z[util])], fmt='o', color='darkred')
    ax[0].fill_between(X_theta, Y_lab_diff - Y_std_lab_diff, Y_lab_diff + Y_std_lab_diff, facecolor='darkred', alpha=0.3)
    ax[0].set_ylabel(r'$\overline{l_{basal} - l_{apical}}$', fontsize=15)
    ax[0].set_xlabel('Position (z/L)', fontsize=15)
    
    ax[1].plot(Y_lab_diff, Y_std_lab_diff, color='darkred')
    ax[1].scatter(lab_diff_mean, std_lab_diff, color='darkred')
    ax[1].set_ylabel(r'$\sigma_{l_{basal} - l_{apical}}$', fontsize=15)
    ax[1].set_xlabel(r'$\overline{l_{basal} - l_{apical}}$', fontsize=15)
    
    der = derivative(np.linspace(z[util].min(), z[util].max(), len(Y_lab_diff)), pars[0], pars[1], pars[3])
    Y_stderror_lab_diff = Y_std_lab_diff/der
    der = derivative(z[util], pars[0], pars[1], pars[3])
    std_error_lab_diff = std_lab_diff[0:len(z[util])]/der
    ax[2].plot(X_theta, Y_stderror_lab_diff, color='darkred')
    ax[2].scatter(z[util], std_error_lab_diff, color='darkred')
    ax[2].set_ylabel(r'$\sigma_{z}(z)$', fontsize=15)
    ax[2].set_xlabel('Position (z/L)', fontsize=15)

    fig.suptitle(f'Length flat side = {l}, σ = {sigma} (over 100 samples)')
    fig.savefig(f'{dir_path}/A0cases/{A0}/DifferentShape/DifferentLength/{l}/sigma/{sigma}/PI-Lab_diff.jpg', bbox_inches='tight')

    lab_ratio_mean = np.mean(lab_ratio, axis=0)
    std_lab_ratio = np.std(lab_ratio, axis=0)

    fig, ax = plt.subplots(1, 3, figsize=(20, 5))
    pars, cov = curve_fit(f=sigmoid, xdata=z[util], ydata=lab_ratio_mean[0:len(z[util])], bounds=(-np.inf, np.inf), maxfev=10**8)
    Y_lab_ratio = sigmoid(np.arange(min(z[util]), max(z[util]), 10e-4), *pars)
    X_theta = np.linspace(z[util].min(), z[util].max(), len(Y_lab_ratio))
    X_Y_cubic_lab_diff = interp1d(z[util], std_lab_ratio[0:len(z[util])], kind='linear')
    Y_std_lab_rario = X_Y_cubic_lab_diff(X_theta)
    ax[0].plot(np.arange(min(z[util]), max(z[util]), 10e-4), Y_lab_ratio, color='forestgreen')
    ax[0].errorbar(z[util], lab_ratio_mean[0:len(z[util])], yerr=std_lab_ratio[0:len(z[util])], fmt='o', color='forestgreen')
    ax[0].fill_between(X_theta, Y_lab_ratio - Y_std_lab_rario, Y_lab_ratio + Y_std_lab_rario, facecolor='forestgreen', alpha=0.3)
    ax[0].set_ylabel(r'$\overline{\frac{l_{apical}}{l_{basal}}}$', fontsize=15)
    ax[0].set_xlabel('Position (z/L)', fontsize=15)
    
    ax[1].plot(Y_lab_ratio, Y_std_lab_rario, color='forestgreen')
    ax[1].scatter(lab_ratio_mean, std_lab_ratio, color='forestgreen')
    ax[1].set_ylabel(r'$\sigma_{\frac{l_{apical}}{l_{basal}}}$', fontsize=15)
    ax[1].set_xlabel(r'$\overline{\frac{l_{apical}}{l_{basal}}}$', fontsize=15)
    
    der = derivative(np.linspace(z[util].min(), z[util].max(), len(Y_std_lab_rario)), pars[0], pars[1], pars[3])
    Y_stderror_lab_ratio = Y_std_lab_rario/der
    der = derivative(z[util], pars[0], pars[1], pars[3])
    std_error_lab_ratio = std_lab_ratio[0:len(z[util])]/der
    ax[2].plot(X_theta, Y_stderror_lab_ratio, color='forestgreen')
    ax[2].scatter(z[util], std_error_lab_ratio, color='forestgreen')
    ax[2].set_ylabel(r'$\sigma_{z}(z)$', fontsize=15)
    ax[2].set_xlabel('Position (z/L)', fontsize=15)
    
    fig.suptitle(f'Length flat side = {l}, σ = {sigma} (over 100 samples)')
    fig.savefig(f'{dir_path}/A0cases/{A0}/DifferentShape/DifferentLength/{l}/sigma/{sigma}/PI-Lab_ratio.jpg', bbox_inches='tight')
    
    