from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
from rdkit import DataStructs
from rdkit import Chem
import numpy as np
import os
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

class Observation:
    """
    Allow the computer agent a to make a decision based on the molecule observations
    """
    
    def __init__(self, mol):
        """
        This is the constructor

        :param mol: The current molecule being observed
        :type mol: RWMol
        """
        self.mol = mol
        self.observation = Chem.MolToSmiles(mol)
        self.info = []
    
   
    def getInfo(self):
        """
        Get information about the molecule

        :return: morgen fingerprint bits and feature information both as numpy.Array objects
        :rtype: tuple
        """
        self.info.clear()
        feats = factory.GetFeaturesForMol(self.mol)
        fp = AllChem.GetMorganFingerprintAsBitVect(self.mol,2,nBits=1024)
        fp_arr = np.zeros((1,))
        for y in feats:
            self.info.append(y.GetType())
        
        DataStructs.ConvertToNumpyArray(fp,fp_arr)
        self.bits = np.nonzero(fp_arr)   
        return self.bits,self.info

    
    def getObservation(self):
        """
        Gets the observations from the molecule

        :return: Smiles string of the current molecule
        :rtype: string
        """
        self.observation = Chem.MolToSmiles(self.mol)
        return self.observation
        
    def update(self,mol):
        """
        Updates the molecule being observed

        :param mol: The current molecule being observed
        :type mol: RWMol
        """
        self.mol = mol 