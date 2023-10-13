import numpy as np
import math

######## Backbone Energy evaluation for closed geomeyry like Round Organoid and the geometry in the Lutolf paper ########


class Energy():
    """Args:
        Object Config
        apical_tension (float): line tension of the apical surface
        lateral_tension (float): line tension of the lateral surface
        compress_cell (float): compressibility modulus of a cell
        A0 (float): preferred area of the cell
    """                 
    def __init__(self, **kwargs):
        
        
        if 'N_cells' in kwargs:
            self.N_cells = kwargs['N_cells']
        else:
            self.N_cells = 0

        if 'apical' in kwargs:
            self.apical = kwargs['apical']
            del kwargs['apical']
        else:
            self.apical = np.zeros((1, 2))

        if 'basal' in kwargs:
            self.basal = kwargs['basal']
            del kwargs['basal']
        else:
            self.basal = np.zeros((1, 2))

        if 'apical_tension' in kwargs:
            self.apical_tension = kwargs['apical_tension']
            del kwargs['apical_tension']
        else:
            self.apical_tension = 1.0

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
            self.A0 = np.zeros(self.N_cells)
        
        if 'dt' in kwargs:
            self.dt = kwargs['dt']
            del kwargs['dt']
        else:
            self.dt = 0.025
        
        if 'seed' in kwargs:
            self.seed = kwargs['seed']
            del kwargs['seed']
        else:
            self.seed = 0
            
        self.E = 0.0

    def evaluate_areas(self):
        areas = np.zeros((self.N_cells, 1))
        j = 0
        for vertex in range(self.N_cells):
            x_coord = np.array([self.apical[vertex, 0], self.apical[(vertex + 1) % self.N_cells, 0], self.basal[(vertex + 1) % self.N_cells, 0], self.basal[vertex, 0]])
            y_coord = np.array([self.apical[vertex, 1], self.apical[(vertex + 1) % self.N_cells, 1], self.basal[(vertex + 1) % self.N_cells, 1], self.basal[vertex, 1]])
            A_cell = 0
            for i in range(len(x_coord)):
                A_cell += x_coord[i]*y_coord[(i+1) % len(x_coord)] - x_coord[(i+1) % len(x_coord)]*y_coord[i]
            areas[j] = np.abs(A_cell)*0.5
            j+=1
        return areas

    def energy(self):
        self.E = 0
        for vertex in range(self.N_cells):
            if self.apical_tension > 0:
                self.E += self.apical_tension*math.sqrt((self.apical[vertex, 0] - self.apical[(vertex + 1) % self.N_cells, 0])**2 + (self.apical[vertex, 1] - self.apical[(vertex + 1) % self.N_cells, 1])**2)
            
            if self.lateral_tension > 0:
                self.E += self.lateral_tension*math.sqrt((self.apical[vertex, 0] - self.basal[vertex, 0])**2 + (self.apical[vertex, 1] - self.basal[vertex, 1])**2)
            
            if self.compress_cell > 0:
                x_coord = np.array([self.apical[vertex, 0], self.apical[(vertex + 1) % self.N_cells, 0], self.basal[(vertex + 1) % self.N_cells, 0], self.basal[vertex, 0]])
                y_coord = np.array([self.apical[vertex, 1], self.apical[(vertex + 1) % self.N_cells, 1], self.basal[(vertex + 1) % self.N_cells, 1], self.basal[vertex, 1]])
                A_cell = 0
                for i in range(len(x_coord)):
                    A_cell += x_coord[i]*y_coord[(i+1) % len(x_coord)] - x_coord[(i+1) % len(x_coord)]*y_coord[i]
                self.E += self.compress_cell*(np.abs(A_cell)*0.5 - self.A0[vertex])**2
            
    
    def distances(self, proposed_apical_x, proposed_apical_y, proposed_basal_x, proposed_basal_y, vertex, side):
        
        dist_lateral = math.sqrt((proposed_apical_x -  proposed_basal_x)**2 + (proposed_apical_y - proposed_basal_y)**2)
    
        if side == 'apical':
            dist_sx = math.sqrt((proposed_apical_x - self.apical[(vertex - 1) % self.N_cells, 0])**2 + (proposed_apical_y - self.apical[(vertex - 1) % self.N_cells, 1])**2)
            dist_dx = math.sqrt((proposed_apical_x - self.apical[(vertex + 1) % self.N_cells, 0])**2 + (proposed_apical_y - self.apical[(vertex + 1) % self.N_cells, 1])**2)      
        else:
            dist_sx = math.sqrt((proposed_basal_x - self.basal[(vertex - 1) % self.N_cells, 0])**2 + (proposed_basal_y - self.basal[(vertex - 1) % self.N_cells, 1])**2)
            dist_dx = math.sqrt((proposed_basal_x - self.basal[(vertex + 1) % self.N_cells, 0])**2 + (proposed_basal_y - self.basal[(vertex + 1) % self.N_cells, 1])**2)
    
        return dist_sx, dist_dx, dist_lateral

    def forces(self, side, vertex):
        fx = np.zeros((self.N_cells, 1))
        fy =  np.zeros((self.N_cells, 1))
        dist_sx, dist_dx, dist_lateral = self.distances(self.apical[vertex, 0], self.apical[vertex, 1], self.basal[vertex, 0], self.basal[vertex, 1], vertex, side)
        if side == 'apical':
                     
            fx[vertex] = -self.lateral_tension*((self.apical[vertex, 0] - self.basal[vertex, 0]) / dist_lateral)
            fy[vertex] = -self.lateral_tension*((self.apical[vertex, 1] - self.basal[vertex, 1]) / dist_lateral)

            fx[vertex] -= self.apical_tension*((self.apical[vertex, 0] - self.apical[(vertex + 1) % self.N_cells, 0]) / dist_dx)
            fx[vertex] -= self.apical_tension*((self.apical[vertex, 0] - self.apical[(vertex - 1) % self.N_cells, 0]) / dist_sx)
            fy[vertex] -= self.apical_tension*((self.apical[vertex, 1] - self.apical[(vertex + 1) % self.N_cells, 1]) /dist_dx)
            fy[vertex] -= self.apical_tension*((self.apical[vertex, 1] - self.apical[(vertex - 1) % self.N_cells, 1])/ dist_sx)
            
            
            if self.compress_cell > 0:
                x_coord1 = np.array([self.apical[vertex, 0], self.apical[(vertex + 1) % self.N_cells, 0], self.basal[(vertex + 1) % self.N_cells, 0], self.basal[vertex, 0]])
                y_coord1 = np.array([self.apical[vertex, 1], self.apical[(vertex + 1) % self.N_cells, 1], self.basal[(vertex + 1) % self.N_cells, 1], self.basal[vertex, 1]])
                x_coord2 = np.array([self.apical[vertex, 0], self.basal[vertex, 0], self.basal[(vertex - 1) % self.N_cells, 0], self.apical[(vertex - 1) % self.N_cells, 0]])
                y_coord2 = np.array([self.apical[vertex, 1], self.basal[vertex, 1], self.basal[(vertex - 1) % self.N_cells, 1], self.apical[(vertex - 1) % self.N_cells, 1]])
                A_cell2 = 0
                A_cell1 = 0
                for i in range(len(x_coord1)):
                    A_cell1 += x_coord1[i]*y_coord1[(i+1) % len(x_coord1)] - x_coord1[(i+1) % len(x_coord1)]*y_coord1[i]
                    A_cell2 += x_coord2[i]*y_coord2[(i+1) % len(x_coord1)] - x_coord2[(i+1) % len(x_coord1)]*y_coord2[i]
                
                fx[vertex] -= self.compress_cell*(0.5*np.abs(A_cell1) - self.A0[vertex])*(-self.basal[vertex, 1] + self.apical[(vertex+1)%self.N_cells, 1])
                fy[vertex] -= self.compress_cell*(0.5*np.abs(A_cell1) - self.A0[vertex])*(self.basal[vertex, 0] - self.apical[(vertex+1)%self.N_cells, 0])
                fx[vertex] -= self.compress_cell*(0.5*np.abs(A_cell2) - self.A0[(vertex - 1) % self.N_cells])*(self.basal[vertex, 1] - self.apical[(vertex-1)%self.N_cells, 1])
                fy[vertex] -= self.compress_cell*(0.5*np.abs(A_cell2) - self.A0[(vertex - 1) % self.N_cells])*(-self.basal[vertex, 0] + self.apical[(vertex-1)%self.N_cells, 0])

        else:
            fx[vertex] = -self.lateral_tension*((self.basal[vertex, 0] - self.apical[vertex, 0]) / dist_lateral)
            fy[vertex] = -self.lateral_tension*((self.basal[vertex, 1] - self.apical[vertex, 1]) / dist_lateral)
            
            if self.compress_cell > 0:
                x_coord1 = np.array([self.apical[vertex, 0], self.apical[(vertex + 1) % self.N_cells, 0], self.basal[(vertex + 1) % self.N_cells, 0], self.basal[vertex, 0]])
                y_coord1 = np.array([self.apical[vertex, 1], self.apical[(vertex + 1) % self.N_cells, 1], self.basal[(vertex + 1) % self.N_cells, 1], self.basal[vertex, 1]])
                x_coord2 = np.array([self.apical[vertex, 0], self.basal[vertex, 0], self.basal[(vertex - 1) % self.N_cells, 0], self.apical[(vertex - 1) % self.N_cells, 0]])
                y_coord2 = np.array([self.apical[vertex, 1], self.basal[vertex, 1], self.basal[(vertex - 1) % self.N_cells, 1], self.apical[(vertex - 1) % self.N_cells, 1]])
                A_cell2 = 0
                A_cell1 = 0
                for i in range(len(x_coord1)):
                    A_cell1 += x_coord1[i]*y_coord1[(i+1) % len(x_coord1)] - x_coord1[(i+1) % len(x_coord1)]*y_coord1[i]
                    A_cell2 += x_coord2[i]*y_coord2[(i+1) % len(x_coord1)] - x_coord2[(i+1) % len(x_coord1)]*y_coord2[i]
                
                fx[vertex] -= self.compress_cell*(0.5*np.abs(A_cell1) - self.A0[vertex])*(-self.basal[(vertex+1)%self.N_cells, 1] + self.apical[vertex, 1])
                fy[vertex] -= self.compress_cell*(0.5*np.abs(A_cell1) - self.A0[vertex])*(self.basal[(vertex+1)%self.N_cells, 0] - self.apical[vertex, 0])
                fx[vertex] -= self.compress_cell*(0.5*np.abs(A_cell2) - self.A0[(vertex - 1) % self.N_cells])*(-self.apical[vertex, 1] + self.basal[(vertex-1)%self.N_cells, 1])
                fy[vertex] -= self.compress_cell*(0.5*np.abs(A_cell2) - self.A0[(vertex - 1) % self.N_cells])*(self.apical[vertex, 0] - self.basal[(vertex-1)%self.N_cells, 0])
                
        return fx, fy
    
    def advance(self):
        mov_x_apical = np.zeros((self.N_cells, 1))
        mov_y_apical = np.zeros((self.N_cells, 1))
        mov_x_basal = np.zeros((self.N_cells, 1))
        mov_y_basal = np.zeros((self.N_cells, 1))
        for vertex in range(self.N_cells):
            mov_x_apical[vertex] = float(self.apical[vertex, 0])
            mov_y_apical[vertex] = float(self.apical[vertex, 1])
            mov_x_basal[vertex] = float(self.basal[vertex, 0])
            mov_y_basal[vertex] = float(self.basal[vertex, 1])

        for vertex in np.random.permutation(self.N_cells):
            #apical
            fx, fy = self.forces('apical', vertex)
            mov_apical_x = float(self.apical[vertex, 0] + self.dt*fx[vertex] / self.motility + self.sigma*np.random.normal(loc=0.0, scale=np.sqrt(self.dt)))
            mov_apical_y = float(self.apical[vertex, 1] + self.dt*fy[vertex] / self.motility + self.sigma*np.random.normal(loc=0.0, scale=np.sqrt(self.dt)))

            mov_x_apical[vertex] = mov_apical_x
            mov_y_apical[vertex] = mov_apical_y
            
            #basal
            fx, fy = self.forces('basal', vertex)
            mov_x= float(self.basal[vertex, 0] + self.dt*fx[vertex] / self.motility + self.sigma*np.random.normal(loc=0.0, scale=np.sqrt(self.dt)))
            mov_y = float(self.basal[vertex, 1] + self.dt*fy[vertex] / self.motility +  self.sigma*np.random.normal(loc=0.0, scale=np.sqrt(self.dt)))

            if mov_y >= 0:
                mov_basal_y = self.r_cup*np.sin(np.arctan(np.abs(mov_y / mov_x)))
                if mov_x < 0:
                    mov_basal_x = -self.r_cup*np.cos(np.arctan(np.abs(mov_y / mov_x)))
                else:
                    mov_basal_x = self.r_cup*np.cos(np.arctan(mov_y / mov_x))
            elif mov_y <= -self.L:
                mov_basal_y = -self.r_cup*np.sin(np.arctan(np.abs((mov_y+self.L)  / mov_x))) - self.L
                if mov_x < 0:
                    mov_basal_x = -self.r_cup*np.cos(np.arctan(np.abs((mov_y+self.L)  / mov_x)))
                else:
                    mov_basal_x = self.r_cup*np.cos(np.arctan(np.abs((mov_y + self.L) / mov_x)))
            else:
                r = math.sqrt(mov_x**2 + mov_y**2)
                mov_basal_x = self.basal[vertex, 0]
                mov_basal_y = -r*np.sin(np.arctan(np.abs(mov_y/mov_x)))
            mov_x_basal[vertex] = mov_basal_x
            mov_y_basal[vertex] = mov_basal_y
            
        self.apical[:, 0] = mov_x_apical[:, 0]
        self.apical[:, 1] = mov_y_apical[:, 0]
        self.basal[:, 0] = mov_x_basal[:,0]
        self.basal[:, 1] = mov_y_basal[:,0]
        self.energy()
        
    def export_dict(self):  
        """Export dictionary containing system's parameters"""
        return {
            "Number of cells": self.N_cells,
            "Motility": self.motility,
            "dt":self.dt
        }