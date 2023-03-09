import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from designspacediscovery import utils as ut
import numpy as np

def test_fingerprint_suite():
    smiles_list = ['CC(C)CC1=CC=C(C=C1)C(C)C(=O)O']
    n_bits = 2048
    mol = Chem.MolFromSmiles(smiles_list[0])
    correct_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, n_bits)

    correct_fp_arr = np.zeros((1,), int)
    DataStructs.ConvertToNumpyArray(correct_fp, correct_fp_arr)

    fp_list = ut.fp_list_from_smiles_list(smiles_list, n_bits = n_bits)

    assert np.array_equal(fp_list[0], correct_fp_arr), 'Mismatch in fingerprints'