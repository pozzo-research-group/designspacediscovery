from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Fragments
import time
import numpy as np
import itertools


def reject_path(paths):
    """
    Find the path that includes the two others

    This function is to be used with the conjugated bridging screening pipeline's bridge 
    finding algorithm. This method naively picks one atom in each ring system, and finds 
    paths between each atom. Unless the molecule has some macro-ring thing, one of these paths
    will include the other two bridges. This function currently works by rejecting the 
    longest path if there are more than 1 paths. A more sophisticated method could find
    the superset of the other ones 

    Parameters:
    -----------
    paths - list of tuples, where each tuple contains the atom indices of a inter-ring bridge path

    Returns:
    --------
    paths - list of tuples, with duplicate path deleted 
    """
    workpath = paths.copy()
    if len(paths) > 1:
        lens = [len(path) for path in workpath]
        del workpath[np.argmax(lens)]
    return workpath

def screen_allowable_atoms(mol, allowable_set = {6, 7, 8, 15, 16, 1, 5}):
    """
    Determine if a molecule is composed of only atoms from a given set. 

    Parameters:
    ----------
    mol - rdkit molecule
    allowable set (set) - set of atomic numbers that are OK

    Returns:
    -------
    bool
    """
    atomset = set([atom.GetAtomicNum() for atom in mol.GetAtoms()])
    return atomset.issubset(allowable_set)

def isRingAromatic(mol, bondRing):
        """
        Determine if a ring is aromatic

        Works by evaluating aromaticity status of every bond in a ring. If all are 
        aromatic, ring is aromatic.

        Parameters:
        ----------
        mol (rdkit mol)
        bondRing (tuple) - Tuple of bond indices for a ring to be evaluated. Get this with RingInfo.BondRings() method.

        Returns:
        -------
        bool
        """
        for ind in bondRing:
            if not mol.GetBondWithIdx(ind).GetIsAromatic():
                return False
        return True

def count_nmember_ring(mol, ri, n = 6, require_aromatic = False):
    """
    Count the number of rings in a molecule with n atoms in them

    Parameters:
    -----------
    mol (rdkit mol)
    ri (rdkit.Chem.rdchem.RingInfo) - ring information object
    n (int) - size of rings to count
    require_aromatic (bool) - whether or not a ring needs to be aromatic to count towards total

    Returns:
    --------
    count (int)
    """
    count = 0
    for ring, bonds in zip(ri.AtomRings(), ri.BondRings()):
        if len(ring) == n:
            if not require_aromatic:
                count+=1
            if require_aromatic:
                if isRingAromatic(mol, bonds):
                    count+=1
    return count


def GetRingSystems(mol, includeSpiro=False):
    """
    Find the independent ring systems in a molecule

    Direct copy-paste out of rdkit cookbook. Finds the independent ring systems, IE chunks of connected rings. Returns that atoms that make up each system as separate tuples

    Parameters:
    ----------
    mol - rdkit mol
    includeSpiro (bool) - whether Spiro systems count as 1 ring system or 2

    Returns:
    -------
    nSystems - list of sets, where each set element is an atom index in a ring system, and each set is one ring system. 
    
    """
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon and (includeSpiro or nInCommon>1):
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    return systems


def check_conjugated_bridging(mol, ring_systems, threshold = 1):
    """
    For molecules with multiple ring systems, checks if connection between systems follows upconverting projects critieria for conjugation.

    Molecule should have more than 1 ring system including spiros. This works by:
    1. Picking an atom in each ring system
    2. Finding paths between each pair of atoms, which should include every possible bridge between rings
    3. Finding the bonds that are not in a ring for those paths
    4. Checking the conjugation status of those bonds
    5. Evaluating our user-defined critieria: Bridging atom is Boron considered conjugated, 

    Parameters:
    ----------
    mol - rdkit molecule
    ring_systems - The ring systems as returned from GetRingSystems. Should be run with spiro = true
    threshold (float) - fraction of bonds in bridges required to be conjugated to return True

    Returns:
    -------
    status (bool) - whether the molecule 'passes' or not
    """
    assert mol is not None, 'Mol should not be None'
    
    n_systems = len(ring_systems)

    key_atoms = []
    
    # randomly pick one atom from each ring system to serve as the endpoint for the path search
    for systems in ring_systems:
        key_atoms.append(np.random.choice(list(systems)))

    paths = []
    # get the shortest path between each pair of atoms
    for atom_1, atom_2 in itertools.combinations(key_atoms, 2):
        try:
            path = Chem.rdmolops.GetShortestPath(mol, int(atom_1), int(atom_2))
        except RuntimeError:
            return False
        paths.append(path)
        #except:
            #fail molecule
    # if there are no paths, probably have some sort of salt with separate rings, don't want to get into that
    if len(paths) == 0:
        return False
    elif any([len(path)==0 for path in paths]):
        return False 
    

    true_paths = reject_path(paths)

    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()

    bridge_conj_status = [] # track the conjugation status of every non-ring bond
    # loop over each of the connections between ring systems
    for path in true_paths:
        nonloop_atoms = [] # track the atoms that aren't in a ring
        for index in path:
            atom = atoms[index]
            if not atom.IsInRing():
                nonloop_atoms.append(atom)
        
        # if there are no non-loop atoms (ie a single carbon bridging two rings
        nonloop_bonds_idx = set() # track the bonds that are not part of a ring. Assume each bond to an atom not in a ring is a non-ring bond
        conjugation_status = [] # track which bonds are not conjugated out of the non-loop set
        if len(nonloop_atoms) == 0:
            # iterate over every combination of ring systems
            for sysA, sysB in itertools.combinations(ring_systems, 2):
                sysA_bonds = set()
                sysB_bonds = set()
                # get system A bonds
                for atomidx in sysA:
                    Abonds = atoms[atomidx].GetBonds()
                    for bond in Abonds:
                        sysA_bonds.add(bond.GetIdx())
                # get system B bonds
                for atomidx in sysB:
                    Bbonds = atoms[atomidx].GetBonds()
                    for bond in Bbonds:
                        sysB_bonds.add(bond.GetIdx())

                # the intersect of these two sets of bond IDs should be the shared bond
                shared_bonds = sysA_bonds.intersection(sysB_bonds)
                nonloop_bonds_idx.update(shared_bonds)
            
        
        # this only works if there is a non-loop atom
        else:
            for atom in nonloop_atoms:
                atombonds = atom.GetBonds()
                # if a nonring atom is boron or phosphorous, just add its bonds to the conjugated set due to our definitions
                if atom.GetAtomicNum() in [5, 15]:
                    for bond in atombonds:
                        conjugation_status.append(True)

                # if not boron proceed as usual
                else:
                    for bond in atombonds:
                        nonloop_bonds_idx.add(bond.GetIdx())

        
        for idx in nonloop_bonds_idx:
            conjugation_status.append(bonds[idx].GetIsConjugated())

        bridge_conj_status.extend(conjugation_status)


    conjcount  = 0
    for entry in bridge_conj_status:
        if entry:
            conjcount +=1

    if conjcount/len(bridge_conj_status) >= threshold:
        status = True
    else:
        status = False
    

    return status
        
def screen_molecule(entry):
    """
    Implements the screening criteria for this project into one function call.

    Expects a dictionary with a 'MolecularWeight' and 'CanonicalSMILES' key.

    Returns a boolean
    """
    
    try:
        molwt = entry['MolecularWeight']
        molwt = float(molwt)
    except KeyError:
        return False
    
    if molwt > 300:
        return False
    try:
        smiles = entry['CanonicalSMILES']
    except KeyError:
        return False
        
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        return False
    
    if mol is None:
        return False
        
    
    #carboxy
    n_carboxyl = Chem.Fragments.fr_C_O(mol)
    n_sulfyrlthingies = Chem.Fragments.fr_C_S(mol)
    if n_carboxyl > 0 or n_sulfyrlthingies > 0:
        return False
    elif mol.HasSubstructMatch(Chem.MolFromSmarts('[cX3]=[OX1]')):
        return False

    
    # atom membership
    if not screen_allowable_atoms(mol):
        return False
    
    #conjugated bridging
    ring_sys = GetRingSystems(mol, includeSpiro=True)
    if len(ring_sys) > 1:
        conjbridge = check_conjugated_bridging(mol, ring_sys)
        if not conjbridge:
            return False
    else:
        pass
    
    # number of aromatic rings
    ri = mol.GetRingInfo()
    count_5 = count_nmember_ring(mol, ri, n = 5, require_aromatic=True)
    count_6 = count_nmember_ring(mol, ri, n = 6, require_aromatic = True)
    
    if (count_5 + count_6) < 2:
        return False
    
    return True