from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import designspacediscovery.upconverting as up

def test_reject_path():
    # 1. test that a set with one path gets returned unmodified
    paths_1 = [(1,2,3)]

    returned_paths = up.reject_path(paths_1)
    assert paths_1[0] == returned_paths[0], '1 path rejection case failed'
    #2. test that longest path successfully removed for >1 path in set
    paths_2 = [(1,2,3), (1,2,3,4,5,6)]
    returned_paths = up.reject_path(paths_2)

    assert returned_paths == [(1,2,3)], 'up.reject_path not rejecting paths correctly'
    return


def test_screen_allowable_atoms():

    #1. test that something with only allowable atoms is passed
    smiles_ok1 = 'CC1=CC(=C(C=C1)PC2=C(C=CC(=C2)OC)O)CO'
    smiles_ok2 = 'CC(C)C1=CSC(=N1)C2=NN=C(O2)C3=CC=CC=C3'

    mol_ok1 = Chem.MolFromSmiles(smiles_ok1)
    mol_ok2 = Chem.MolFromSmiles(smiles_ok2)

    assert up.screen_allowable_atoms(mol_ok1), 'OK atoms rejected'
    assert up.screen_allowable_atoms(mol_ok2), 'OK atoms rejected'

    #2. Test that something with F fails
    smiles_fail = 'CC1(C(C1C(=O)OC(C#N)C2=CC(=C(C=C2)F)OC3=CC=CC=C3)C=C(Cl)Cl)C'

    mol_fail = Chem.MolFromSmiles(smiles_fail)
    assert not up.screen_allowable_atoms(mol_fail), 'Bad molecule passed atom screen'
    return

def test_isRingAromatic():

    #1. Test that an aromatic ring passes
    aromatic_smiles = "C1=CC=CC=C1"
    aromatic_mol = Chem.MolFromSmiles(aromatic_smiles)
    rings = aromatic_mol.GetRingInfo()
    bondring = rings.BondRings()[0]
    assert up.isRingAromatic(aromatic_mol, bondring), 'up.isRingAromatic failed on benzene'

    #2. test that a non-aromatic ring returns false
    #cyclohexane
    aryl_smiles = 'C1CCCCC1'
    aryl_mol = Chem.MolFromSmiles(aryl_smiles)
    rings = aryl_mol.GetRingInfo()
    bondring = rings.BondRings()[0]
    assert not up.isRingAromatic(aryl_mol, bondring), 'up.isRingAromatic thinks cyclohexane aromatic'

    return

def test_count_nmember_ring():

    #1. Test 5 member ring, aromatic
    smiles_5 = 'C1CCCC1'
    mol_5 = Chem.MolFromSmiles(smiles_5)
    ring_5 = mol_5.GetRingInfo()

    assert up.count_nmember_ring(mol_5, ring_5, n = 5, require_aromatic=False) == 1, 'Miscount on cyclopentane rings'



    #2. test 6 memeber ring, nonaromatic
    smiles_6 = 'C1CCCCC1'
    mol_6 = Chem.MolFromSmiles(smiles_6)
    ring_6 = mol_6.GetRingInfo()

    assert up.count_nmember_ring(mol_6, ring_6, n = 6, require_aromatic = False) == 1, 'Miscount on cyclopentane rings'
    #3. test 6 member ring, either aromatic or not
    #anthracene
    smiles_aromatic = 'C1=CC=C2C=C3C=CC=CC3=CC2=C1'
    mol_aromatic = Chem.MolFromSmiles(smiles_aromatic)
    ring_aromatic = mol_aromatic.GetRingInfo()

    assert up.count_nmember_ring(mol_aromatic, ring_aromatic, n = 6, require_aromatic= True) == 3, 'Anthracene miscount'

    return


def test_GetRingSystems():

    #1. test something with no rings
    smiles = 'CCCCC'
    assert up.GetRingSystems(Chem.MolFromSmiles(smiles)) == [], 'Found rings in pentane'
    #2. test something with 1 system
    smiles = 'C1CCCCC1'
    mol = Chem.MolFromSmiles(smiles)

    assert up.GetRingSystems(mol) == [{0,1,2,3,4,5}], 'Found wrong ring in cyclohexane'

    #3. test something with 2 systems
    smiles = 'C1=CC=C(C=C1)NC2=CC=CC=C2' #diphenyl amine
    mol = Chem.MolFromSmiles(smiles)

    assert up.GetRingSystems(mol) == [{0, 1, 2, 3, 4, 5}, {7, 8, 9, 10, 11, 12}], 'wrong rings on DPA'

    #4. Test spiro case
    spirosmiles = 'C1CCC2(CC1)CCCCC2'
    spiromol = Chem.MolFromSmiles(spirosmiles)

    assert up.GetRingSystems(spiromol, includeSpiro=False) == [{0, 1, 2, 3, 4, 5}, {3, 6, 7, 8, 9, 10}], 'issue on spiro true case for spiro'
    assert up.GetRingSystems(spiromol, includeSpiro= True) == [{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}], 'Issue on no spiro case on spiro for ring system ID'

    return

def test_check_conjugated_bridging():

    #1. check for dpa
    smiles = 'C1=CC=C(C=C1)NC2=CC=CC=C2' #diphenyl amine
    mol = Chem.MolFromSmiles(smiles)
    ring_systems = up.GetRingSystems(mol, includeSpiro=True)
    assert up.check_conjugated_bridging(mol, ring_systems), 'DPA failed conjugation check'

    #2. Check for dpa but with boron
    smiles = 'C1=CC=C(C=C1)BC2=CC=CC=C2' #diphenyl amine
    mol = Chem.MolFromSmiles(smiles)
    ring_systems = up.GetRingSystems(mol, includeSpiro=True)
    assert up.check_conjugated_bridging(mol, ring_systems), 'Boron bridge failed conjugation check'
    #3. check for something with a non-conjugated bridge
    smiles = 'C1=CC=C(C=C1)CCCC2=CC=CC=C2'
    mol = Chem.MolFromSmiles(smiles)
    ring_systems = up.GetRingSystems(mol, includeSpiro=True)
    assert not up.check_conjugated_bridging(mol, ring_systems), 'aryl bridge passed conjugation check but should not have'

    # try a conjugated carbon bridge

    smiles = 'C1=CC=C(C=C1)C=CC2=CC=CC=C2' #conjugated carbon bridge
    mol = Chem.MolFromSmiles(smiles)    
    ring_systems = up.GetRingSystems(mol, includeSpiro=True)
    assert up.check_conjugated_bridging(mol, ring_systems), 'conjugated bridge failed conjugation check'
    return




