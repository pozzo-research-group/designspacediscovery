# functions related to running 2d similarity search
import designspacediscovery.querypubchem as qpc

def is_integery(val):
    """
    Figure out if something could be a pubchem id
    """
    if (not isinstance(val, str)) and (not isinstance(val, int)):
        return False
    else:
        if isinstance(val, str):
            try:
                int(val)
                return True
            except:
                return False
        elif isinstance(val, int):
            return True
        else:
            return False



def find_similar_molecules(basis_set:dict, threshold = 90, max_records = 5000, representation = 'CID')-> dict:
    """
    Find similar molecules using the pubchem fast2dsimilarity api

    Uses tanimoto similarity scores based on pubchem fingerprints

    Parameters:
    -----------
    basis_set (dict): dictionary of molecules with desired property in format {key:CID}. Molecules must be represented as pubchem CID or SMILES

    """ 
    assert isinstance(basis_set, dict), 'basis set must be a dictionary with Pubchem CIDs as values'
    assert is_integery(list(basis_set.values())[0]), 'Basis set values must be Pubchem CIDs' 

    url_dict = {}
    for key in list(basis_set.keys()):
        cid = basis_set[key]
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/{cid}/cids/JSON?Threshold={threshold}&MaxRecords={max_records}'
        url_dict[key] = url
    
    retriever = qpc.pubchemQuery()
    similarity_responses = retriever.run_queries(url_dict)

    similarities  = {k:v.json()['IdentifierList']['CID'] for k, v in similarity_responses.items() if v != 'FAILED'}

    return similarities
