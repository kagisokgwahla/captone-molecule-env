import gym
from gym import error, spaces, utils
from gym.utils import seeding
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import RWMol 
from rdkit import RDConfig
from rdkit import DataStructs
import numpy as np
import os
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

ADD     = "add"
REMOVE  = "remove"
FRONT   = "front"
BACK    = "back"