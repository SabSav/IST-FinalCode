import numpy as np
import math

class Config:
    """Args:
        N_cells (int): total number of cell
        h (float): height of the external monolayer
        apical (array): 1D array (tot_cells, ) of positions. 
        basal (array): 1D array of positions. 
    """
    def __init__(self, N_cells=10, alpha=10, beta=8, h=5, wavelength = 0.25, **kwargs):

        self.N_cells= N_cells
        self.alpha = alpha
        self.h = h
        self.beta = beta
        self.wavelength = wavelength
        assert self.alpha + self.h > self.beta

        if 'basal' in kwargs:
            self.basal = np.asarray(kwargs['basal'])
            del kwargs['basal']
        else:
            self.basal = np.zeros((self.N_cells+1, 2))
            self.basal[0, 0] = 0
            for i in range(1, self.N_cells):
                self.basal[i, 0] = i
            self.basal[self.N_cells, 0] = self.N_cells
            self.basal[:, 1] = self.beta*np.cos(2*np.pi*self.basal[:, 0] / (self.N_cells*self.wavelength))
    
        if 'apical' in kwargs:
            self.apical = np.asarray(kwargs['apical'])
            del kwargs['apical']
        else:
            self.apical = np.zeros((self.N_cells+1, 2))
            self.apical[0, 0] = 0
            for i in range(1, self.N_cells):
                self.apical[i, 0] = i
            self.apical[self.N_cells, 0] = self.N_cells
            self.apical[:, 1] = self.alpha*np.cos(2*np.pi*self.apical[:, 0] / (self.N_cells*self.wavelength)) + self.h
            
        if 'density' in kwargs:
            self.density = kwargs['density']
            del kwargs['density']
        else:
            self.density = self.wavelength*self.N_cells+1

        if 'areas' in kwargs:
            self.areas = kwargs['areas']
            del kwargs['areas']
        else:
            self.areas = np.zeros((self.N_cells, 1))
        
        if 'min_distance' in kwargs:
            self.min_distance = kwargs['min_distance']
            del kwargs
        else:
            self.min_distance = 0.5
        if 'seed' in kwargs:
            self.seed = kwargs['seed']
            del kwargs['seed']
        else:
            self.seed = 0

        if 'Ncells_wavelength' in kwargs:
            self.N_cells_wavelength = kwargs['Ncells_wavelength']
        else:
            self.N_cells_wavelength = 6        

class Energy(Config):
    """Args:
        Object Config
        apical_tension (float): line tension of the apical surface
        lateral_tension (float): line tension of the lateral surface
        compress_cell (float): compressibility modulus of a cell
        A0 (float): preferred area of the cell
    """                 
    def __init__(self, N_cells=10, alpha=10, beta=1.5, **kwargs):
        super().__init__(N_cells, alpha, beta,**kwargs)

        if 'apical_tension' in kwargs:
            self.apical_tension = kwargs['apical_tension']
            del kwargs['apical_tension']
        else:
            self.apical_tension = 1.0
        
        if 'basal_tension' in kwargs:
            self.basal_tension = kwargs['basal_tension']
            del kwargs['basal_tension']
        else:
            self.basal_tension = 1.0

        if 'lateral_tension' in kwargs:
            self.lateral_tension = kwargs['lateral_tension']
            del kwargs['lateral_tension']
        else:
            self.lateral_tension = 1.0

        if 'compress_cell' in kwargs:
            self.compress_cell = kwargs['compress_cell']
            del kwargs['compress_cell']
        else:
            self.compress_cell = 1.0
        
        if 'motility' in kwargs:
            self.motility = kwargs['motility']
            del kwargs['motility']
        else:
            self.motility = 1.0
        
        if 'A0' in kwargs:
            self.A0 = kwargs['A0']
            del kwargs['A0']
        else:
            self.A0 = np.zeros((self.N_cells, 1))
        
        if 'dt' in kwargs:
            self.dt = kwargs['dt']
            del kwargs['dt']
        else:
            self.dt = 0.025

        self.E = 0.0

    def evaluate_areas(self):
        areas = np.zeros((self.N_cells, 1))
        j = 0
        for vertex in range(self.N_cells):
            x_coord = np.array([self.basal[vertex, 0], self.basal[vertex+1, 0], self.apical[vertex+1, 0], self.apical[vertex, 0]])
            y_coord = np.array([self.basal[vertex, 1], self.basal[vertex+1, 1], self.apical[vertex+1, 1], self.apical[vertex, 1]])
            A_cell = 0
            for i in range(len(x_coord)):
                A_cell += x_coord[i]*y_coord[(i+1) % len(x_coord)] - x_coord[(i+1) % len(x_coord)]*y_coord[i]
            areas[j] = np.abs(A_cell)*0.5
            j+=1

        return areas
    
    def energy(self):
        self.E = 0
        for vertex in range(self.N_cells+1): #in a not closed layer Ncells means Ncells + 1 verteces
            #apical
            if self.apical_tension > 0:
                if vertex < self.N_cells: #to avoid double counting when arriving at the last cell
                    self.E += self.apical_tension*math.sqrt((self.apical[vertex, 0] - self.apical[(vertex + 1), 0])**2 + (self.apical[vertex, 1] - self.apical[(vertex + 1), 1])**2)
            
            if self.lateral_tension > 0:
                self.E += self.lateral_tension*math.sqrt((self.apical[vertex, 0] - self.basal[vertex, 0])**2 + (self.apical[vertex, 1] - self.basal[vertex, 1])**2)
            
            if self.compress_cell > 0:
                if vertex < self.N_cells:
                    x_coord = np.array([self.basal[vertex, 0], self.basal[vertex+1, 0], self.apical[vertex+1, 0], self.apical[vertex, 0]])
                    y_coord = np.array([self.basal[vertex, 1], self.basal[vertex+1, 1], self.apical[vertex+1, 1], self.apical[vertex, 1]])
                    A_cell = 0
                    for i in range(len(x_coord)):
                        A_cell += x_coord[i]*y_coord[(i+1) % len(x_coord)] - x_coord[(i+1) % len(x_coord)]*y_coord[i]
                    self.E += self.compress_cell*(np.abs(A_cell)*0.5 - self.A0[vertex])**2


    def distances(self, proposed_apical_x, proposed_apical_y, proposed_basal_x, proposed_basal_y, vertex, side):
        
        dist_lateral = math.sqrt((proposed_apical_x -  proposed_basal_x)**2 + (proposed_apical_y - proposed_basal_y)**2)
    
        if side == 'apical':
            dist_sx = math.sqrt((proposed_apical_x - self.apical[(vertex - 1), 0])**2 + (proposed_apical_y - self.apical[(vertex - 1), 1])**2)
            dist_dx = math.sqrt((proposed_apical_x - self.apical[(vertex + 1), 0])**2 + (self.apical[(vertex + 1), 1] - proposed_apical_y)**2)
            
        else:
            dist_sx = math.sqrt((proposed_basal_x - self.basal[(vertex - 1), 0])**2 + (proposed_basal_y - self.basal[(vertex - 1), 1])**2)
            dist_dx = math.sqrt((proposed_basal_x - self.basal[(vertex + 1), 0])**2 + (self.basal[(vertex + 1), 1] - proposed_basal_y)**2)

        return dist_sx, dist_dx, dist_lateral

    def forces(self, side, vertex):
        fx = np.zeros((self.N_cells+1, 1))
        fy =  np.zeros((self.N_cells+1, 1))
        dist_sx, dist_dx, dist_lateral = self.distances(self.apical[vertex, 0], self.apical[vertex, 1], self.basal[vertex, 0], self.basal[vertex, 1], vertex, side)
        if side == 'apical':            
            fx[vertex] = -self.lateral_tension*((self.apical[vertex, 0] - self.basal[vertex, 0]) / dist_lateral)
            fy[vertex] = -self.lateral_tension*((self.apical[vertex, 1] - self.basal[vertex, 1]) / dist_lateral)
        
            fx[vertex] -= self.apical_tension*((self.apical[vertex, 0] - self.apical[(vertex + 1), 0]) / dist_dx)
            fx[vertex] -= self.apical_tension*((self.apical[vertex, 0] - self.apical[(vertex - 1), 0]) / dist_sx)
            fy[vertex] -= self.apical_tension*((self.apical[vertex, 1] - self.apical[(vertex + 1), 1]) /dist_dx)
            fy[vertex] -= self.apical_tension*((self.apical[vertex, 1] - self.apical[(vertex - 1), 1])/ dist_sx)
            
            if self.compress_cell > 0:
                x_coord1 = np.array([self.basal[vertex, 0], self.basal[vertex+1, 0], self.apical[vertex+1, 0], self.apical[vertex, 0]])
                y_coord1 = np.array([self.basal[vertex, 1], self.basal[vertex+1, 1], self.apical[vertex+1, 1], self.apical[vertex, 1]])
                x_coord2 = np.array([self.basal[vertex, 0], self.apical[vertex, 0], self.apical[vertex-1, 0], self.basal[vertex-1, 0]])
                y_coord2 = np.array([self.basal[vertex, 1], self.apical[vertex, 1], self.apical[vertex-1, 1], self.basal[vertex-1, 1]])
                A_cell2 = 0
                A_cell1 = 0
                for i in range(len(x_coord1)):
                    A_cell1 += x_coord1[i]*y_coord1[(i+1) % len(x_coord1)] - x_coord1[(i+1) % len(x_coord1)]*y_coord1[i]
                    A_cell2 += x_coord2[i]*y_coord2[(i+1) % len(x_coord1)] - x_coord2[(i+1) % len(x_coord1)]*y_coord2[i]
            
                fx[vertex] -= self.compress_cell*(0.5*np.abs(A_cell1) - self.A0[vertex])*(self.basal[vertex, 1] - self.apical[vertex+1, 1])
                fy[vertex] -= self.compress_cell*(0.5*np.abs(A_cell1) - self.A0[vertex])*(-self.basal[vertex, 0] + self.apical[vertex+1, 0])
                fx[vertex] -= self.compress_cell*(0.5*np.abs(A_cell2) - self.A0[vertex-1])*(-self.basal[vertex, 1] + self.apical[vertex-1, 1])
                fy[vertex] -= self.compress_cell*(0.5*np.abs(A_cell2) - self.A0[vertex-1])*(self.basal[vertex, 0] - self.apical[vertex-1, 0])
            
        else:
            fx[vertex] = -self.lateral_tension*((self.basal[vertex, 0] - self.apical[vertex, 0]) / dist_lateral)
        
            if self.compress_cell > 0:
                x_coord1 = np.array([self.basal[vertex, 0], self.basal[vertex+1, 0], self.apical[vertex+1, 0], self.apical[vertex, 0]])
                y_coord1 = np.array([self.basal[vertex, 1], self.basal[vertex+1, 1], self.apical[vertex+1, 1], self.apical[vertex, 1]])
                x_coord2 = np.array([self.basal[vertex, 0], self.apical[vertex, 0], self.apical[vertex-1, 0], self.basal[vertex-1, 0]])
                y_coord2 = np.array([self.basal[vertex, 1], self.apical[vertex, 1], self.apical[vertex-1, 1], self.basal[vertex-1, 1]])
                A_cell2 = 0
                A_cell1 = 0
                for i in range(len(x_coord1)):
                    A_cell1 += x_coord1[i]*y_coord1[(i+1) % len(x_coord1)] - x_coord1[(i+1) % len(x_coord1)]*y_coord1[i]
                    A_cell2 += x_coord2[i]*y_coord2[(i+1) % len(x_coord1)] - x_coord2[(i+1) % len(x_coord1)]*y_coord2[i]
                
                fx[vertex] -= self.compress_cell*(0.5*np.abs(A_cell1) - self.A0[vertex])*(self.basal[vertex+1, 1] - self.apical[vertex, 1])
                fx[vertex] -= self.compress_cell*(0.5*np.abs(A_cell2) - self.A0[vertex-1])*(-self.basal[vertex-1, 1] + self.apical[vertex, 1])
            
        return fx, fy

        
    def advance(self):
        mov_x_apical = np.zeros((self.N_cells+1, 1))
        mov_y_apical = np.zeros((self.N_cells+1, 1))
        mov_x_basal = np.zeros((self.N_cells+1, 1))
        mov_y_basal = np.zeros((self.N_cells+1, 1))
        #for vertex in np.random.permutation(np.arange(2, self.N_cells-1)): #confluency condition: first and last cell are fixed
        for vertex in np.random.permutation(np.arange(1, self.N_cells)): #when there is noise the permutation might be useful for avoiding correlation between sequential vertices movements (?)
            #apical
            fx, fy = self.forces('apical', vertex)
            for i in range(0, self.N_cells_wavelength):
                if i*self.wavelength*self.N_cells <= self.apical[vertex, 0] <= (i + 1)*self.wavelength*self.N_cells: #density = Ncells/wavelength kept constant
                    if i*self.wavelength*self.N_cells <= self.apical[vertex, 0] + self.dt*fx[vertex] / self.motility <= (i + 1)*self.wavelength*self.N_cells:
                        mov_x_apical[vertex] = self.apical[vertex, 0] + self.dt*fx[vertex] / self.motility 
                        mov_y_apical[vertex] = self.apical[vertex, 1] + self.dt*fy[vertex] / self.motility
                        i = self.N_cells_wavelength
                    else:
                        mov_x_apical[vertex] = self.apical[vertex, 0]
                        mov_y_apical[vertex] = self.apical[vertex, 1]
                        i = self.N_cells_wavelength

            #basal
            fx, fy = self.forces('basal', vertex)
            for i in range(0, self.N_cells_wavelength):
                if i*self.wavelength*self.N_cells <= self.basal[vertex, 0] <= (i + 1)*self.wavelength*self.N_cells:
                    if i*self.wavelength*self.N_cells <= self.basal[vertex, 0] + self.dt*fx[vertex] / self.motility <= (i + 1)*self.wavelength*self.N_cells:
                        mov_x_basal[vertex] = self.basal[vertex, 0] + self.dt*fx[vertex] / self.motility
                        mov_y_basal[vertex] =  self.beta*np.cos(2*np.pi*mov_x_basal[vertex] / (self.N_cells*self.wavelength)) 
                        i = self.N_cells_wavelength
                    else:
                        mov_x_basal[vertex] = self.basal[vertex, 0]
                        mov_y_basal[vertex] = self.basal[vertex, 1]
                        i = self.N_cells_wavelength
        self.apical[1:self.N_cells, 0] = mov_x_apical[1:self.N_cells, 0]
        self.apical[1:self.N_cells, 1] = mov_y_apical[1:self.N_cells, 0]
        self.basal[1:self.N_cells, 0] = mov_x_basal[1:self.N_cells, 0]
        self.basal[1:self.N_cells, 1] = mov_y_basal[1:self.N_cells, 0]
        self.energy()

            

    def export_dict(self):
        """Export dictionary containing system's parameters"""
        return {
            "Number of cells": self.N_cells,
            "Motility": self.motility
        }
            