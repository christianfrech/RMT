import numpy as np
import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import linalg as LA
import math
import xlsxwriter
from matplotlib import cm

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

#Creates the list of upper limits for bonding distances (functions identically to np.linspace)
def createUpperLimitList(low,high,spacing):
    limitlist=[]
    for i in range(int((high-low)/spacing)+1):
        value=low+(spacing*i)
        limitlist.append(value)
    return limitlist

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
        
    def get_atom_elneg(self):
        if os.path.isfile(self.file):
            pdb_file = open(self.file, 'r')
        else:
            raise OSError('File {} does not exist'.format(self.file))
        
        lineno = 0
        frames = []
        atoms = []
        
        for line in pdb_file:
            lineno += 1
            if line.startswith('ATOM'):
                 at_obj = PDBAtom(line)
                 atoms.append(at_obj.electronegativity)
            elif line.startswith('END'):
                    frames.append(atoms)
                    atoms = []
        pdb_file.close()
        
        base = np.zeros((len(frames), len(frames[0])))
        for i in range(len(frames)):
            for j in range(len(frames[i])):
                base[i, j] = frames[i][j]
        self.elneg_adj = np.reshape(base, (len(frames), 1, len(frames[0]))) - np.reshape(base, (len(frames), len(frames[0]), 1))
        
        return self.elneg_adj
    
    def get_elneg_adj(self, t = 4):
        if len(self.elneg_adj) == 1:
            self.get_atom_elneg()
        self.elneg_adj = (self.elneg_adj < t).astype(int)
        return self.elneg_adj
    
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

        self.adjacency_graphs = ((self.distance_graphs < t) & (self.distance_graphs > s)).astype(int)

        #Eliminating same-atom intramolecular interactions:
        for frame in range(len(self.adjacency_graphs)):  
            hydrogenbond_count=0
            for i in range(len(self.adjacency_graphs[frame])):
                for j in range(len(self.adjacency_graphs[frame][i])):
                    if (i//3==j//3):
                        self.adjacency_graphs[frame][i][j]=1
                    elif ((self.valence_list[frame][i]!=self.valence_list[frame][j]) & (self.adjacency_graphs[frame][i][j]==1)):
                        #self.adjacency_graphs[frame][i][j]=1
                        if (i,j) not in used_valence_list:
                            hydrogenbond_count+=1
                        else: 
                            pass
                        used_valence_list.add((i,j))
                        used_valence_list.add((j,i))
                    else:
                        self.adjacency_graphs[frame][i][j]=0

        hydrogenbonds_array.append(hydrogenbond_count)

        return self.adjacency_graphs
    
    def get_elec_adj(self):
        #if len(self.adjacency_graphs) == 1:
        self.get_atom_adj(lowerlimit,upperlimit)
            
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
        return self.elec_adjacency_graphs
    
    def make_eigenvalues(self, hamiltonian_iter=10):
        #self.elec_adjacency_graphs=[]
        self.elec_adjacency_graphs=self.get_elec_adj()
        elec_count = len(self.elec_adjacency_graphs[0])
        self.eigenvalues = np.zeros((len(self.elec_adjacency_graphs), hamiltonian_iter, elec_count))
        for frame in range(len(self.elec_adjacency_graphs)):
            frame_eigs = []
            for i in range(hamiltonian_iter):
                print(i)
                r = np.random.normal(size=(elec_count, elec_count))
                rt = np.transpose(r)
                h = (r + rt) / np.sqrt(2 * elec_count)          
                adj_r = self.elec_adjacency_graphs[frame] * h
                eigs = np.ndarray.tolist(np.linalg.eigvals(adj_r))
                eigs.sort()
                #for i in range(len(eigs)):
                frame_eigs.append(np.real(eigs))
            self.eigenvalues[frame] = frame_eigs
        return self.eigenvalues

    def get_spacings(self,types):
        allframes_spacing=[]
        eigenvalues = self.make_eigenvalues()
        #eigenvalues=[np.random.rand(1000,60)]
        eigenvalues=np.array(eigenvalues)
        if types=='all':
            medians=[]
            nexttomedian=[]
            spacings=[]
            allspacings=np.zeros((len(eigenvalues[0]),len(eigenvalues[0][0])-1))
            for i in range(len(eigenvalues[0])):
                for j in range(len(eigenvalues[0][0])-1):
                    allspacings[i][j]=(eigenvalues[0][i][j+1]-eigenvalues[0][i][j])
            allspacings=allspacings.transpose()
            counts, binedges = np.histogram(allspacings)

            return [counts, allspacings]
        else:
            for frame in range(len(self.elec_adjacency_graphs)):
                medians=[]
                nexttomedian=[] 
                spacings=[]
                for i in range(len(eigenvalues[frame])):
                    #Calculating medians and converting eigenvalue array to 1xn list:
                    medians.append(np.median(eigenvalues[frame][i]))
                    if len(eigenvalues[frame][i])%2==0:
                        nexttomedian.append((eigenvalues[frame][i][math.floor(len(eigenvalues[frame])/2)+1]+eigenvalues[frame][i][math.floor(len(eigenvalues[frame])/2)+2])/2)
                    else:
                        nexttomedian.append(eigenvalues[frame][i][math.floor(len(eigenvalues[frame])/2)+1])
                #Calculating median and median+1 spacings, along with spacing standard deviation:
                for i in range(len(medians)):
                    spacings.append(abs(nexttomedian[i]-medians[i]))
                allframes_spacing.append(spacings)
            return allframes_spacing

    def get_spacings_and_stdevs(self,types,framerange):
        eigenvalues = self.make_eigenvalues()
        spacings=[]
        medians=[]
        nexttomedian=[]
        stdevs=[]
        eigenvalues=np.array(eigenvalues)

        '''
        Choose lower allspacings code block for first five eigenvalue spacings, upper allspacings code block for all eigenvalue spacings
        '''
        if types=='all':
            '''
            allspacings=np.zeros((len(eigenvalues[frame]),len(eigenvalues[frame][0])-1))
            for i in range(len(eigenvalues[frame])):
                for j in range(5):
                    allspacings[i][j]=(eigenvalues[frame][i][j+1]-eigenvalues[frame][i][j])
            '''
            for frame in framerange:
                allspacings=np.zeros((len(eigenvalues[frame]),5))
                for i in range(len(eigenvalues[frame])):
                    for j in len(eigenvalues[frame][i]):
                        allspacings[i][j]=(eigenvalues[frame][i][j+1]-eigenvalues[frame][i][j])

                allspacings=allspacings.transpose()
                for i in range(len(allspacings)):
                    stdevs.append(np.std(allspacings[i]))
            return [stdevs, allspacings]

        else:
            stdevs=[]
            for frame in range(len(self.elec_adjacency_graphs)):
                spacings=[]
                medians=[]
                nexttomedian=[]
                for i in range(len(eigenvalues[frame])):
                #Calculating medians and converting eigenvalue array to 1xn list:
                    medians.append(np.median(eigenvalues[frame][i]))
                    if len(eigenvalues[frame][i])%2==0:
                        nexttomedian.append((eigenvalues[frame][i][math.floor(len(eigenvalues[frame])/2)+1]+eigenvalues[frame][i][math.floor(len(eigenvalues[frame])/2)+2])/2)
                    else:
                        nexttomedian.append(eigenvalues[frame][i][math.floor(len(eigenvalues[frame])/2)+1])
                #Calculating median and median+1 spacings, along with spacing standard deviation:
                for i in range(len(medians)):
                    spacings.append(nexttomedian[i]-medians[i])
                stdevs.append(np.std(spacings))
            return stdevs

    def get_stdevs(self):
        stdevs=[]
        medianvalues=[]
        eigenvalues = self.make_eigenvalues()
        for i in range(len(eigenvalues)):
            eigenvalues[i].sort()
            medianvalues.append(eigenvalues[i])
        stdevs.append(np.std(medianvalues))
        return stdevs


'''
Initialize the frames and distances to iterate through:
lowerlimit/upperlimit: lower and upper limit of bonding distances
frameindices: which frames are used for stdev calculation
types: which type of eigenvalue spacings to use 
    ('all_frames'=3D surface plot of bond distances (X) vs. frame number (Y) vs. standard deviation of spacings (Z))
    ('all'=3D surface plot of bond distances (X) vs. eigenvalue pair number (Y) vs. standard deviation of spacings (Z))
    ('median'=2D line plot of bond distances (X) vs. standard deviation of spacings (Y))
'''
lowerlimit=1
upperlimit_list=createUpperLimitList(lowerlimit,3,1)
allframes_stdevs=[]
framesindices=np.linspace(1,2500,5)
types='all_frames'
frames=[]


'''
Calling methods from Adj_Mats class, iterating through frames (frameindex) and through bonding distances (upperlimit) to
append the standard deviations to full 2D list (allframes_stdevs)
'''

oneframe_stdevs=[]
hydrogenbonds_array=[]
lowerlimit=1

for upperlimit in upperlimit_list:
    if __name__ == "__main__":
        print(upperlimit)
        full = Adj_Mats('water.pdb')
        results=full.get_spacings_and_stdevs(types,framesindices)
        if (types=='all'):
            oneframe_stdevs.append(results[0])
            allspacings=results[1]
        else:
            oneframe_stdevs=results
    else:
        pass

    allframes_stdevs.append(oneframe_stdevs)


allframes_stdevs=np.array(allframes_stdevs)
print(hydrogenbonds_array)

'''
MAKING PLOTS (different plots for different input of types variable)
'''
if (types=='all_frames'):
    Z = allframes_stdevs
    fig = plt.figure(figsize=plt.figaspect(0.5))

    #Creating X-axis for 3D plots (frames)
    xx=framesindices
    iterations=[upperlimit_list,hydrogenbonds_array]
    xaxis_labels=['Cut-Off Distance For Interactions (Angstroms)', 'Number of Hydrogen Bonds']

    for col in range(2):
        if col==0:  
            yy = upperlimit_list
            X, Y = np.meshgrid(xx,yy)
            print(X.shape)
            print(Y.shape)
            print(Z.shape)
            ax = fig.add_subplot(1, 2, 1, projection='3d')
            surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            fig.colorbar(surf, shrink=0.5, aspect=10)
            ax.set_xlabel('Frame index')
            ax.set_ylabel(xaxis_labels[col])
            ax.set_zlabel('Spacings Standard Deviation')
        else:
            yy = hydrogenbonds_array
            X, Y = np.meshgrid(xx,yy)
            ax = fig.add_subplot(1, 2, 2, projection='3d')
            surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            fig.colorbar(surf, shrink=0.5, aspect=10)
            ax.set_xlabel('Frame index')
            ax.set_ylabel(xaxis_labels[col])
            ax.set_zlabel('Spacings Standard Deviation')

    plt.savefig('3DWater_{1}frame_hydrogenbonds_standarddeviation.png', dpi=95)
    plt.show()

elif (types=='all'):
    mins=[]
    maxes=[]
    spacingmin=min(range(5))
    spacingmax=max(range(5))
    spacingiterations=len(range(5))
    bondmin=lowerlimit
    bondmax=max(upperlimit_list)
    bonditerations=len(upperlimit_list)

    eigenvaluearray=np.array(allspacings)

    #XY-plane
    yy=np.linspace(0,bondmax,bonditerations)
    xx=np.linspace(spacingmin,spacingmax,spacingiterations)
    X,Y=np.meshgrid(xx,yy)

    #Z value Distributions
    allstdevs=np.array(oneframe_stdevs)
    Z = np.array(allstdevs)
    print(X.shape)
    print(Y.shape)
    print(Z.shape)

    #Make 3D surface plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z)

    ax.set_xlabel('Eigenvalue Spacing Pair #')
    ax.set_ylabel('Maximum Hydrogen Bond Length (Angstroms)')
    ax.set_zlabel('Standard Deviation of Spacings')

    plt.savefig('3D_watermolecule_standarddeviation_allspacings_hbondlength.png', dpi=95)
    plt.show()

else:
    ###PLOTS###
    list_stdevs=list(allframes_stdevs[0])
    fig, axes = plt.subplots(1,2)

    bondmin = lowerlimit
    bondmax = max(upperlimit_list)
    bondspacings = len(upperlimit_list)
    bonditeration = np.linspace(bondmin, bondmax, bondspacings)

    iterations=[bonditeration,hydrogenbonds_array]
    xaxis_labels=['Cut-Off Distance For Interactions (Angstroms)', 'Number of Hydrogen Bonds']

    print(list_stdevs)


    for col in range(len(axes)):
        axes[col].plot(iterations[col], list_stdevs)  
        axes[col].set_xlabel(xaxis_labels[col])
        axes[col].set_ylabel('Spacings Standard Deviation')


    plt.savefig('3DWater_{1}frame_hydrogenbonds_standarddeviation.png', dpi=95)
    plt.show()


#NEXT STEPS:
#WITHOUT INTERMOLECULAR RESONANCE:
    #Take 2500 frames and make adjacencies of all of them
    #2) Identify only covalent bonds in all of them 
    #3) Identify only hydrogen bonds
    #4) INTRAmolecular interaction are resonating, INTERmolecular interactions only dissimilar bonds (O-H)
    #5) For spacing graph, we want minimum bond case, maximum bond case, and middle bond case (spacing vs probability of spacing) <---- move further right, lower stdev
    #6) (Eigenvalue distributiion vs probability of value)--minimum bond case, maximum bond case, and middle bond case  <---- become more semicircular
    #7) Standard deviation vs. distance for bonds AND THEN # h-bonds  <---- should decrease
    #8) (Entropy vs dist of bonding) AND THEN (Entropy vs # of bonds) <---- should increase


    # DO EVERYTHING ABOVE FOR FIRST AND LAST TIME FRAME

    #Read about Structured Random Matrices