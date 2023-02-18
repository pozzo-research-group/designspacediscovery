# functions related to running 2d similarity search
import designspacediscovery.querypubchem as qpc

def find_similar_molecules(basis_set:dict, threshold = 90, max_records = 5000, representation = 'CID')-> dict:
    """
    Find similar molecules using the pubchem fast2dsimilarity api

    Uses tanimoto similarity scores based on pubchem fingerprints

    Parameters:
    -----------
    basis_set (dict): dictionary of molecules with desired property in format {key:CID}. Molecules must be represented as pubchem CID or SMILES
    


    """ 