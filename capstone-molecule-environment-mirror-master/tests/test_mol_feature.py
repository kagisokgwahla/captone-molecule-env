import unittest
from .context import mol_feature

import numpy as np
class TestMolFeature(unittest.TestCase):

    def setUp(self):
        mol = 'C(=O)OC1=CC=CC=C1'
        self.feat = mol_feature.Mol_Feature(mol)
       
    def test_getSmile(self):
        self.assertEquals(self.feat.getSmile(), 'C(=O)OC1=CC=CC=C1')

    def test_contains(self):
        self.assertTrue(self.feat.contains(np.array(['Arom6'])))
        self.assertFalse(self.feat.contains(np.array(['Arom5'])))

        self.assertTrue(self.feat.contains(np.array(['Arom6','RH6_6'])))
        self.assertFalse(self.feat.contains(np.array(['Arom5','RH5_5 '])))

        self.assertTrue(self.feat.contains(np.array(['SingleAtomDonor' ,'SingleAtomAcceptor', 'SingleAtomAcceptor' ,'AcidicGroup'])))
        self.assertFalse(self.feat.contains(np.array(['Imidazole','ZnBinder1','BasicGroup','Guanidine'])))

if __name__ == '__main__':
    unittest.main()
