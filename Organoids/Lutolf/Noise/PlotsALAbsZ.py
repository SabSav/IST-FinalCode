import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import os
import json
import numpy as np
import math 

if __name__ == '__main__':  
    

    def sigmoid(x, a,b,c,d):
        return c + d*(1/ (1 + np.exp(-x/b + a)))
   
    sigmas = [0.01, 0.05, 0.1, 0.2]
    Nsamples = 100
    dir_path = os.path.dirname(os.path.realpath(__file__))
    l = [7.5, 10, 15, 20]
    Ncells = [20, 24, 30, 36]
    dt = 1e-3
    rcup = 2.5
    for sig in sigmas:
        ds_angles = []
        lmean_sample_apical = []
        var_sample_apical = []
        for leng in range(len(l)):
            file = open(f'{dir_path}/DifferentLength{l[leng]}/sigma/{sig}/bundle-{0}.json')
            data = json.load(file)
            lapical = np.zeros((Nsamples, Ncells[leng]))
            ds_angle = []
            Nf = data['Number of cells on flat'] + 1
            Nc = int(0.5*(0.5*Ncells[leng] - Nf))
            theta = np.linspace(0, 0.5*np.pi, Nc)
            for i in range(len(theta)):
                ds_angle.append(np.cos(theta[i])*rcup)
            vertex = len(ds_angle) - 1
            for i in range(vertex + 1, vertex + 1 + Nf):
                ds_angle.append(ds_angle[i - 1] - l[leng]/Nf)
            theta = np.flip(theta)
            for i in range(1, len(theta)):
                ds_angle.append(-l[leng] - np.cos(theta[i])*rcup)
            ds_angle = np.array([0.5*(i - max(ds_angle))/(min(ds_angle) - max(ds_angle)) for i in ds_angle])
            ds_angles.append(ds_angle)

            for sample in range(Nsamples):
                file = open(f'{dir_path}/{l[leng]}/sigma/{sig}/bundle-{sample}.json')
                data = json.load(file)
                apical_y = data['Apical_y'][-1]
                apical_x = data['Apical_x'][-1]
                basal_x = data['Basal_x'][-1]
                basal_y = data['Basal_y'][-1]
                for vertex in range(Ncells[leng]):
                    lapical[sample, vertex] = math.sqrt((apical_x[vertex] - apical_x[(vertex+1)%Ncells[leng]])**2 + (apical_y[vertex] - apical_y[(vertex+1)%Ncells[leng]])**2)
            lmean_sample_apical.append(np.mean(lapical, axis=0)) #mean for each vertex
            var_sample_apical.append(np.var(lapical, axis=0)) #std for each vertex
        
        lapical_to_plot = []
        std_to_plot = []
        #lapical_to_plot.append(lapical[0])
        for leng in range(len(l)):
            vertex = 0
            lapical_to_store = []
            std_to_store = []
            for pos0 in range(len(ds_angles[0])):
                lapical_aggr = 0
                var_agrr = 0
                if vertex < len(ds_angles[leng]):
                    if ds_angles[leng][vertex] <= ds_angles[0][pos0]:
                        lapical_aggr += lmean_sample_apical[leng][vertex]
                        var_agrr +=  var_sample_apical[leng][vertex]
                        vertex += 1
                    std_to_store.append(np.sqrt(var_agrr/vertex))
                    lapical_to_store.append(lapical_aggr)
            lapical_to_plot.append(lapical_to_store)
            std_to_plot.append(std_to_store)

        lapical_to_plot = np.array(lapical_to_plot)
        std_to_plot = np.array(std_to_plot)
        vertex_util_cr = []
        vertex_util_dcr = []
        for vertex in range(len(ds_angles[0])):
            if 0 <= ds_angles[0][vertex] <= 0.25:
                vertex_util_cr.append(vertex)
            else:
                vertex_util_dcr.append(vertex)
        
        color = []
        np.random.seed(0)
        for i in range(len(l)):
            color.append('#%06X' % np.random.randint(0, 0xFFFFFF))
        plt.figure(figsize=(8,8))
        for i in range(len(l)):
            X_theta = np.linspace(ds_angles[0][vertex_util_cr].min(), ds_angles[0][vertex_util_cr].max(), 500)
            X_Y_cubic_lapical = interp1d(ds_angles[0][vertex_util_cr], lapical_to_plot[i][vertex_util_cr], kind='linear')
            X_Y_cubic_std = interp1d(ds_angles[0][vertex_util_cr], std_to_plot[i][vertex_util_cr], kind='linear')
            Y_lapical = X_Y_cubic_lapical(X_theta)
            Y_std = X_Y_cubic_std(X_theta)
            plt.errorbar(ds_angles[0][vertex_util_cr]*2*(2*rcup + l[0]), lapical_to_plot[i][vertex_util_cr], yerr=std_to_plot[i][vertex_util_cr], color=color[i], fmt='o')    
            plt.fill_between(X_theta*2*(2*rcup + l[0]), Y_lapical - Y_std, Y_lapical + Y_std, facecolor=color[i], alpha=0.2)
            pars, cov = curve_fit(f=sigmoid, xdata=(ds_angles[0][vertex_util_cr])*2*(2*rcup + l[0]), ydata=lapical_to_plot[i][vertex_util_cr], bounds=(-np.inf, np.inf), maxfev=10**8)
            plt.plot(np.arange(min((ds_angles[0][vertex_util_cr])*(2*rcup + l[0])), max((ds_angles[0][vertex_util_cr])*2*(2*rcup + l[0])), 10e-4), sigmoid(np.arange(min((ds_angles[0][vertex_util_cr])*2*(2*rcup + l[0])), max((ds_angles[0][vertex_util_cr])*2*(2*rcup + l[0])), 10e-4), *pars),
            color=color[i], linestyle='solid', linewidth=2, label=f'L = {l[i]}. a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f' % tuple(pars))
            
            file_parameters = {'Sigma': sig}
            file_parameters['a'] = pars[0]
            file_parameters['b'] = pars[1]
            file_parameters['c'] = pars[2]
            file_parameters['d'] = pars[3]
            output = f'{dir_path}/{l[i]}/sigma/{sig}/file_parameters_angle.json'
            with open(output, 'w') as file:
                json.dump(file_parameters, file, indent=20)
        plt.legend()
        plt.title(f'Noise σ = {sig}')
        plt.xlabel('z', fontsize=15)
        plt.ylabel(r'$l_{apical}(z)$', fontsize=15)
        plt.ylim(0, 3.0)
        plt.savefig(f'{dir_path}/ApicalLengthFunction-{sig}.jpg')