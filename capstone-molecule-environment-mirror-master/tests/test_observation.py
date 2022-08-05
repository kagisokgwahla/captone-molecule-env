import unittest
from .context import observation
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import RWMol
from rdkit import RDConfig
import numpy as np
import os
fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

class TestObservation(unittest.TestCase):

    def setUp(self):
        self.mol = RWMol(Chem.MolFromSmiles("C"))
        self.obs = observation.Observation(self.mol)

    def tearDown(self):
        pass
    
    def test_getInfo(self):
        info = []
        bits = []
        feats = factory.GetFeaturesForMol(self.mol)
        fp = AllChem.GetMorganFingerprintAsBitVect(self.mol,2,nBits = 1024)
         
        for y in feats:
            info.append(y.GetType())
        fp_arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp,fp_arr)
        bits = np.nonzero(fp_arr)   
       

        self.assertEqual(self.obs.getInfo(), (bits,info))

    def test_getObservation(self):
        observations = Chem.MolToSmiles(self.obs.mol)

        self.assertEquals(self.obs.getObservation(), observations)

    def test_update(self):
        mol = RWMol(Chem.MolFromSmiles("CC"))
        self.obs.update(mol)

        self.assertEquals(self.obs.mol, mol)


if __name__ == '__main__':
    unittest.main()
