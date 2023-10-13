from Energy import Energy
import numpy as np
import math

######## Advance step function for the different closed geometries ########

class advanceRound(Energy):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if 'sigma' in kwargs:
            self.sigma = kwargs['sigma']
        else:
            self.sigma = 0

        if 'r' in kwargs:
            self.r = kwargs['r']
        else:
            self.r = 0
        
        if 'h' in kwargs:
            self.h = kwargs['r']
        else:
            self.h = 0

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

        for vertex in np.random.permutation(self.N_cells): #when there is noise the permutation might be useful for avoiding correlation between sequential vertices movements (?)
            #apical vertices: free to move
            fx, fy = self.forces('apical', vertex)

            mov_x_apical[vertex] = float(self.apical[vertex, 0] + self.dt*fx[vertex] / self.motility + self.sigma*np.random.normal(loc=0.0, scale=np.sqrt(self.dt)))
            mov_y_apical[vertex] = float(self.apical[vertex, 1] + self.dt*fy[vertex] / self.motility + self.sigma*np.random.normal(loc=0.0, scale=np.sqrt(self.dt)))
            
            #basal vertices: need to be constrained on the inital circle, it is just a proojection
            fx, fy = self.forces('basal', vertex)
            mov_x= float(self.basal[vertex, 0] + self.dt*fx[vertex] / self.motility + self.sigma*np.random.normal(loc=0.0, scale=np.sqrt(self.dt)))
            mov_y = float(self.basal[vertex, 1] + self.dt*fy[vertex] / self.motility +  self.sigma*np.random.normal(loc=0.0, scale=np.sqrt(self.dt)))

            if mov_x >= 0:
                mov_basal_x = self.r*np.cos(np.arctan(np.abs(mov_y/mov_x)))
            else:
                mov_basal_x = -self.r*np.cos(np.arctan(np.abs(mov_y/mov_x)))
            if mov_y >= 0:
                mov_basal_y = self.r*np.sin(np.arctan(np.abs(mov_y/mov_x)))
            else:
                mov_basal_y = -self.r*np.sin(np.arctan(np.abs(mov_y/mov_x)))

            mov_x_basal[vertex] = mov_basal_x
            mov_y_basal[vertex] = mov_basal_y
            
        self.apical[:, 0] = mov_x_apical[:, 0]
        self.apical[:, 1] = mov_y_apical[:, 0]
        self.basal[:, 0] = mov_x_basal[:,0]
        self.basal[:, 1] = mov_y_basal[:,0]
        self.energy()

class advanceLutolfExperiment(Energy):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if 'r_cup' in kwargs:
            self.r_cup = kwargs['r_cup']
        else:
            self.r_cup = 0.25
        
        if 'N_cells_flat' in kwargs:
            self.N_cells_flat = kwargs['N_cells_flat']
        else:
            self.N_cells_flat = 0

        if 'h' in kwargs:
            self.h = kwargs['h']
        else:
            self.h = 0.5
        
        if 'L' in kwargs:
            self.L = kwargs['L']
        else:
            self.L = 5.0

        if 'sigma' in kwargs:
            self.sigma = kwargs['sigma']
        else:
            self.sigma = 0
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

        for vertex in np.random.permutation(self.N_cells): #when there is noise the permutation might be useful for avoiding correlation between sequential vertices movements (?)
            #apical vertices: free to move
            fx, fy = self.forces('apical', vertex)
            mov_apical_x = float(self.apical[vertex, 0] + self.dt*fx[vertex] / self.motility + self.sigma*np.random.normal(loc=0.0, scale=np.sqrt(self.dt)))
            mov_apical_y = float(self.apical[vertex, 1] + self.dt*fy[vertex] / self.motility + self.sigma*np.random.normal(loc=0.0, scale=np.sqrt(self.dt)))

            mov_x_apical[vertex] = mov_apical_x
            mov_y_apical[vertex] = mov_apical_y
            
            #basal vertices: need to be constrained on the inital shape, it is just a projection
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