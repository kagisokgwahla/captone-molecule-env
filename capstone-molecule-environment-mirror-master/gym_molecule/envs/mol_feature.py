import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import RWMol 
from rdkit import RDConfig
from rdkit import DataStructs
import os
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)


class Mol_Feature:
    """
    Provides easy access to molecule features
    """

    def __init__(self,smiles):
        """
        This is the constructor

        :param smiles: The smile string for the molecule
        :type smiles: string
        """

        self.smiles = smiles
        self.mol = RWMol(Chem.MolFromSmiles(smiles))
        
        #create a feature a numpy array
        self.feats = factory.GetFeaturesForMol(self.mol)
        self.feature_arr = np.array([y.GetType() for y in self.feats])
        print(self.feature_arr)
        
        #create a morgen finger print array 
        self.fp = AllChem.GetMorganFingerprintAsBitVect(self.mol,2,nBits=1024)
        self.fp_arr = np.zeros((1,0))
        DataStructs.ConvertToNumpyArray(self.fp,self.fp_arr)
        np.nonzero(self.fp_arr)
        
    def getSmile(self):
        """
        Used to get the smiles string of the molecule

        :return: The smiles string for the molecule
        :rtype: string
        """
        return self.smiles
    
    def contains(self,query):
        """
        Used to query molecule features and morgen fingerprints

        :param query: An array molecule features and morgen fingerprint to be compared
        :type query: numpy.Array
        :return: Whether atleast one of the item in the query array is present in the molecules features or morgen fingerprints
        :rtype: bool
        """
        # check if query contains value in feature array  print list
        if (len(self.feature_arr) !=0) & (query.size != 0) :
            for feature in self.feature_arr:
                for item in query:
                    if feature == item:
                        return True
                
        # check if query contains value in morgen fingerprint array
        if (self.fp_arr.size !=0) & (query.size != 0) :
            for fp in self.fp_arr:
                for item in query:
                    if item == fp:
                        return True
        
                    
        return False