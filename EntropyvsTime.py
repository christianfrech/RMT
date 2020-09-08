import numpy as np
import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import linalg as LA
import xlsxwriter
import math

# Fills transpose of mat[N][N] in tr[N][N] 
def transpose(mat, tr, N): 
    for i in range(N): 
        for j in range(N): 
            tr[i][j] = mat[j][i] 
   
# Returns true if mat[N][N] is symmetric, else false 
def isSymmetric(mat, N): 
      
    tr = [ [0 for j in range(len(mat[0])) ] for i in range(len(mat)) ] 
    transpose(mat, tr, N) 
    for i in range(N): 
        for j in range(N): 
            if (mat[i][j] != tr[i][j]): 
                return False
    return True

#Make a matrix symmetric
def makeSymmetric(inputmatrix):
    inputmatrixT = np.transpose(inputmatrix)
    inputmatrix = (inputmatrix + inputmatrixT)
    '''
    for i in range(len(inputmatrix)):
        for j in range(len(inputmatrix[0])):
            if inputmatrix[i][j] not in [0,1]:
                inputmatrix[i][j]=1
    '''
    return inputmatrix

def createUpperLimitList(low,high,spacing):
    limitlist=[]
    for i in range(int((high-low)/spacing)+1):
        value=low+(spacing*i)
        limitlist.append(value)
    return limitlist

def binSizes(values,binscount):
    binvalues=[]
    if framecount==1:
        width=abs(max(values)-min(values))
        for i in range(binscount):
            print(i)
            binvalues.append(sum([min(values),i*(width/binscount)]))
    else:
        for frame in range(len(values)):
            width=abs(max(values[frame])-min(values[frame]))
            for i in range(binscount):
                print(i)
                binvalues.append(sum([min(values[frame]),i*(width/binscount)]))
    return binvalues

class PDBAtom(object):
    """ Class to represent a single atom's position and state at a frame
    
    Attributes:
        _valence_dict (dict{str: int}): A dictionary of valence electron count per element
        x (float): The x coordinate of the atom
        y (float): The y coordinate of the atom
        z (float): The z coordinate of the atom
        valence_count (int): Number of valence electrons in the atom
    """
    
    
    _valence_dict = {'C': 4,
                     'H': 1,
                     'N': 5,
                     'O': 6,
                     'S': 6}
    
    _electroneg_dict = {'C': 2.55,
                        'H': 2.2,
                        'N': 3.04,
                        'O': 3.44,
                        'S': 2.58}
    
    def __init__(self, string):
        """ Standard PDB file format
        ATOM    277  O1  LYS A  14      21.138  -0.865  -4.761  1.00  0.00           O1-
        """
#       Coordinate Parser 
        self.x = float(string[30:38].strip())
        self.y = float(string[38:46].strip())
        self.z = float(string[46:54].strip())
        
#       Element and Valence Electron Number Parser
        self.element_spec = string[77:].strip()
        mod = 0
        if self.element_spec.endswith(('-', '+')):
            self.element_sym = self.element_spec[:-2].strip()
            mod = int(self.element_spec[-2])
            mod *= (-1, 1)[self.element_spec.endswith('-')]
        else:
            self.element_sym = self.element_spec.strip()
        self.valence_count = PDBAtom._valence_dict.get(self.element_sym)
        if self.valence_count is None:
            raise TypeError('Used an element that is not in the valence dictionary')
        else:
            self.valence_count += mod
        self.electronegativity = PDBAtom._electroneg_dict.get(self.element_sym)


class Adj_Mats(object):
    """ Class to represent a series of adjacency matrices
    
    Attributes:
        file (str): The path of the pdb file to be parsed
        valence_list (Array[Array[int]]): Stores the number of valence electrons in all atoms in every frame
        distance_graphs (Array[Array[Array[int]]]): The series of distance matrices of the atoms in the evolution
        adjacenecy_graphs (Array[Array[Array[int]]]): The series of adjacency matrices of the atoms in the evolution
        elec_adjacency_graphs (Array[Array[Array[int]]]): The series of adjacency matrices of electrons in the evolution
        
    Methods:
        set_atom_dists: Used to set the distance_graphs attribute
        set_atom_adj: Used to set the adjacency_graphs attribute
        get_atom_dists: Used to parse the pdb file to create a distance_graphs object
        get_atom_adj: Used to set an adjacency threshold on the distance matrices and make adjacency matrices
    """
    
    
    def __init__(self, pdb):
        self.file = pdb
        self.valence_list = np.zeros(1, int)
        self.distance_graphs = np.zeros(1, int)
        self.adjacency_graphs = np.zeros(1, int)
        self.elec_adjacency_graphs = np.zeros(1, int)
        self.elneg_adj = np.zeros(1, int)
        self.eigenvalues = None
        self.bin_probs = None
        self.entropy = None
        self.energy = None
        self.cont_ent = None
    
    def set_atom_dists(self, new_dists):
        self.distance_graphs = new_dists
        
    def set_atom_adj(self, new_adj):
        self.adjacency_graphs = new_adj
           
    def get_atom_dists(self):
        if os.path.isfile(self.file):
            pdb_file = open(self.file,'r')
        else:
            raise OSError('File {} does not exist'.format(self.file))

        lineno = 0
        frames = []
        atoms = []
        val_frames = []
        val_atoms = []
        
        for line in pdb_file:
            lineno += 1
            if line.startswith('ATOM'):
                try:
                    at_obj = PDBAtom(line)
                    atoms.append([at_obj.x, at_obj.y, at_obj.z])
                    val_atoms.append(at_obj.valence_count)
                except:
                    sys.stderr.write('\nProblem parsing line {} in file {}\n'.format(lineno, self.file))
                    sys.stderr.write(line)
                    sys.stderr.write('Probably ATOM entry is formatted incorrectly?\n')
                    sys.stderr.write('Please refer to - http://www.wwpdb.org/documentation/format32/sect9.html#ATOM\n\n')
                    sys.exit(1)
            elif line.startswith('END'):
                frames.append(atoms)
                atoms = []
                val_frames.append(val_atoms)
                val_atoms = []
        pdb_file.close()
    
        base = np.zeros((len(framesindices), len(frames[0]), 3))
        for i in range(len(framesindices)):
            for j in range(len(frames[i])):
                for k in range(len(frames[i][j])):
                    base[i][j][k] = frames[i][j][k]
        dists = np.reshape(base, (len(framesindices), 1, len(frames[0]), 3)) - np.reshape(base, (len(framesindices), len(frames[0]), 1, 3))  
        dists = dists**2
        dists = dists.sum(3)
        dists = np.sqrt(dists)
        
        self.valence_list = val_frames
        self.distance_graphs = dists
        
        return self.distance_graphs
    
    def get_atom_adj(self, s=1, t=4):
        if len(self.distance_graphs) == 1:
            self.get_atom_dists()

        used_valence_list = set()
        hydrogenbond_count=0

        self.adjacency_graphs = ((self.distance_graphs < t) & (self.distance_graphs > s)).astype(int)

        #Eliminating same-atom intramolecular interactions:
        for frame in range(len(self.adjacency_graphs)):
            for i in range(len(self.adjacency_graphs[frame])):
                for j in range(len(self.adjacency_graphs[frame][i])):
                    if (i//3==j//3):
                        self.adjacency_graphs[frame][i][j]=1
                    elif ((self.valence_list[frame][i]==self.valence_list[frame][j]) & (i//3!=j//3)):
                        self.adjacency_graphs[frame][i][j]=0
                    elif ((self.valence_list[frame][i]!=self.valence_list[frame][j]) & (self.adjacency_graphs[frame][i][j]==1)):
                        self.adjacency_graphs[frame][i][j]=1
                        if (i,j) not in used_valence_list:
                            hydrogenbond_count+=1
                        else: 
                            pass
                    else:
                        self.adjacency_graphs[frame][i][j]=0

                    used_valence_list.add((i,j))
                    used_valence_list.add((j,i))

        hydrogenbonds_array.append(hydrogenbond_count)

        return self.adjacency_graphs
    
    def get_elec_adj(self):
        #if len(self.adjacency_graphs) == 1:
        self.get_atom_adj()
            
        total_val = 0
        
        for i in range(len(self.valence_list[0])):
            total_val += self.valence_list[0][i]
        valencelistframes = len(self.valence_list)
        self.elec_adjacency_graphs = np.zeros((len(framesindices), total_val, total_val))
        curr_n, curr_m = 0, 0
        
        for i in range(len(self.adjacency_graphs)):
            for j in range(len(self.adjacency_graphs[0])):
                for b in range(self.valence_list[i][j]):
                    for k in range(len(self.adjacency_graphs[0][0])):
                        for a in range(self.valence_list[i][k]):
                            self.elec_adjacency_graphs[i][curr_n][curr_m] = self.adjacency_graphs[i][j][k]
                            curr_m += 1
                    curr_m = 0
                    curr_n += 1
            curr_n = 0

        if (isSymmetric(self.elec_adjacency_graphs[0], len(self.elec_adjacency_graphs[0])))==True: 
            pass
        else: 
            print("No")

        return self.elec_adjacency_graphs
    
    def entropyCalculation(self):
        mat1=self.get_elec_adj()
        #n=number of rows
        n=len(mat1[0])
        #m=off-diagonal entries
        m=0
        batch_entropies=[]

        for frame in range(len(mat1)):
            for i in range(len(mat1[frame])):
                for j in range(len(mat1[frame])):
                    if i!=j:
                        if mat1[0][i][j]!=0:
                            m+=1
                        else: pass
                    else: pass
            entropy = m/2+n/2-(m/2)*math.log(2)+(m/2+n/2)*math.log(n)-(m/2+n/2)*math.log(math.pi)
            batch_entropies.append(entropy)
        #EigenEntropies.append(entropy)
        return batch_entropies


binscount=1000
lowerlimit=1
frame_list=[np.linspace(0,499,500),np.linspace(500,999,500),np.linspace(1000,1499,500),np.linspace(1500,1999,500),np.linspace(2000,499,500)]
hydrogenbonds_array=[]
entropies=[]
hydrogenbonds_oneframe=[]
allframes_entropies=np.array(len(frame_list)*500)

data_folder = "/home/gemsec-user/Desktop/"
file_to_open = data_folder + "water.pdb"
file = open(file_to_open)


for frameindex in frame_list:
    if __name__ == "__main__":
        full = Adj_Mats(file_to_open)
        batch = full.entropyCalculation()
        for i in range(len(allframes_entropies)):
            allframes_entropies[((frameindex*500) + i)] = batch[i]
    else: pass


#naming protocol
names_allframes = []
for j in framesindices:
    name_string = 'Frame ' + str(int(j))
    names_allframes.append(name_string)

print(names_allframes)

###REMOVE [0] FROM PLOT WHEN USING MORE THAN 1 FRAME####

###PLOTS###


entropies_array=np.array(allframes_entropies)
fig, ax = plt.subplots(1,1)

framemin = min(framesindices)
framemax = max(framesindices)
framespacings = len(framesindices)
bonditeration = np.linspace(framemin, framemax, framespacings)

print(entropies_array)
print(hydrogenbonds_array)



ax.plot(bonditeration, entropies_array)
ax.set_xlabel('Frame (10ns spacing)')
ax.set_ylabel('Entropy (nats)')


plt.savefig('3DWater_timesteps_entropyplot.png', dpi=95)
plt.show()