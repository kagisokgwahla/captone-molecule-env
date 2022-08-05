import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import sqrt
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
from rdkit import Chem
from rdkit.Chem import Descriptors


class Datacapture:
    """Processes the training data to compare our current molecule with the drugs in the market
    """
    
    def __init__(self,current_mol):
        """This is the constructor

        :param current_mol: The current molecule
        :type current_mol: RWMol
        """
        self.current_mol = current_mol
        self.sol = pd.read_csv('gym_molecule/envs/resources/delaney.csv')
        self.mol_list= []
        for element in self.sol.SMILES:
              mol = Chem.MolFromSmiles(element)
              self.mol_list.append(mol)
            
   
    def generate(self,smiles, verbose=False):
        """Get required information about all the molecules in the data set

        :param smiles: List of smiles strings of all the molecules in the data set
        :type smiles: string[]
        :param verbose: [description], defaults to False
        :type verbose: bool, optional
        :return: tabular data of all the molecules in the data set
        :rtype: pandas.DataFrame
        """

        moldata= []
        for elem in smiles:
            mol=Chem.MolFromSmiles(elem) 
            moldata.append(mol)
            
        baseData= np.arange(1,1)
        i=0  
        for mol in moldata:        

            desc_MolLogP = Descriptors.MolLogP(mol)
            desc_MolWt = Descriptors.MolWt(mol)
            desc_NumRotatableBonds = Descriptors.NumRotatableBonds(mol)

            row = np.array([desc_MolLogP,
                            desc_MolWt,
                            desc_NumRotatableBonds])   

            if(i==0):
                baseData=row
            else:
                baseData=np.vstack([baseData, row])
            i=i+1      

        columnNames=["MolLogP","MolWt","NumRotatableBonds"]   
        descriptors = pd.DataFrame(data=baseData,columns=columnNames)

        return descriptors


    def current_generate(self):
        """Gets the required information about the current molecule

        :return: Tabular data of the current molecule
        :rtype: pandas.DataFrame
        """

        mol=self.current_mol 

        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_MolWt = Descriptors.MolWt(mol)
        desc_NumRotatableBonds = Descriptors.NumRotatableBonds(mol)

        row = np.array([desc_MolLogP,
                        desc_MolWt,
                        desc_NumRotatableBonds])   

        baseData= row
        baseData=np.vstack([baseData])
            

        columnNames=["MolLogP","MolWt","NumRotatableBonds"]   
        descriptors = pd.DataFrame(data=baseData,columns=columnNames)

        return descriptors
    
    def AromaticAtoms(self,m):
        """Gets the number of aromatic atoms in the molecule

        :param m: The current molecule in the environment
        :type m: RWMol
        :return: The number of aromatic atoms
        :rtype: int
        """
        aromatic_atoms = [m.GetAtomWithIdx(i).GetIsAromatic() for i in range(m.GetNumAtoms())]
        aa_count = []
        for i in aromatic_atoms:
            if i==True:
                aa_count.append(1)
        sum_aa_count = sum(aa_count)
        return sum_aa_count
        
    def processing(self):
        """Plots the graph, trains and evaluates the model
        """
        ######################################## Data preparation #########################################

        df = self.generate(self.sol.SMILES)
        desc_AromaticAtoms = [self.AromaticAtoms(element) for element in self.mol_list]
        desc_HeavyAtomCount = [Descriptors.HeavyAtomCount(element) for element in self.mol_list]
        desc_AromaticProportion = [self.AromaticAtoms(element)/Descriptors.HeavyAtomCount(element) for element in self.mol_list]
        df_desc_AromaticProportion = pd.DataFrame(desc_AromaticProportion, columns=['AromaticProportion'])
        
        curr_df = self.current_generate()
        curr_X = curr_df[['MolLogP', 'MolWt']].values.reshape(-1,2)

        curr_Y = 0.26 - 0.74*(curr_df.MolLogP[0])- 0.0066*(curr_df.MolWt[0])+ 0.0032*(curr_df.NumRotatableBonds[0])-0.4*(self.AromaticAtoms(self.current_mol)/Descriptors.HeavyAtomCount(self.current_mol))

        df = self.generate(self.sol.SMILES)

        X = df[['MolLogP', 'MolWt']].values.reshape(-1,2)
        Y = self.sol.iloc[:,2]

        ######################## Prepare model data point for visualization ###############################
        curr_x = curr_X[:, 0]
        curr_y = curr_X[:, 1]
        curr_z = curr_Y
   

        x = X[:, 0]
        y = X[:, 1]
        z = Y

        
        x_pred = np.linspace(min(x), max(x))   # range of porosity values
        y_pred = np.linspace(min(y), max(y))  # range of brittleness values
        xx_pred, yy_pred = np.meshgrid(x_pred, y_pred)
        model_viz = np.array([xx_pred.flatten(), yy_pred.flatten()]).T


        ################################################ Train #############################################

        ols = linear_model.LinearRegression()
        model = ols.fit(X, Y)
        predicted = model.predict(model_viz)

        ############################################## Evaluate ############################################

        r2 = model.score(X, Y)

        ############################################## Plot ################################################
        
        # calculating centriod
        x_centriod = 0
        y_centriod = 0
        z_centriod = 0
        nx = len(x)
        for xi in x:
            x_centriod += xi
        ny = len(y)
        for yi in y:
            y_centriod += yi    
        nz = len(z)
        for zi in z:
            z_centriod += zi  
            
        x_centriod = x_centriod/nx
        y_centriod = y_centriod/ny
        z_centriod = z_centriod/nz
      

        plt.style.use('default')

        fig = plt.figure()

        ax = fig.add_subplot( projection='3d')
#         ax2 = fig.add_subplot(132, projection='3d')
#         ax3 = fig.add_subplot(133, projection='3d')

#         axes = [ax1, ax2, ax3]
       
        ax.plot(x, y, z, color='k', zorder=15, linestyle='none', marker='o', alpha=0.5)
        ax.scatter(xx_pred.flatten(), yy_pred.flatten(), predicted, facecolor=(0,0,0,0), s=20, edgecolor='#70b3f0')
        ax.set_xlabel('MolLogP', fontsize=12)
        ax.set_ylabel('MolWt', fontsize=12)
        ax.set_zlabel('Predicted LogS', fontsize=12)
        ax.locator_params(nbins=4, axis='x')
        ax.locator_params(nbins=5, axis='x')

        ax.plot(curr_x, curr_y,curr_z,color='r', zorder=15, linestyle='none', marker='o', alpha=1)
        ax.plot([x_centriod], [y_centriod],[z_centriod],color='b', zorder=15, linestyle='none', marker='o', alpha=1)
    
        
        
        # distance from centroid 
        dist = sqrt((curr_x)**2 + (curr_y)**2 + (curr_z)**2)
        
        ax.view_init(elev=60, azim=165)
#         ax2.view_init(elev=4, azim=114)
#         ax3.view_init(elev=60, azim=165)
        
       
        fig.suptitle('$Dist_{mol-centroid} = %.2f$ $  R^2 = %.2f$' % (dist,r2) )
        

        for ii in np.arange(0, 360, 1):
            ax.view_init(elev=32, azim=ii)
            fig.savefig('gym_molecule/envs/resources/%d.png' % ii)
        
