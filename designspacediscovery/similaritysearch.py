# functions related to running 2d similarity search
import designspacediscovery.querypubchem as qpc
import designspacediscovery.utils as utils




def find_similar_molecules(basis_set:dict, threshold = 90, max_records = 5000, representation = 'CID')-> dict:
    """
    Find similar molecules using the pubchem fast2dsimilarity api

    Uses tanimoto similarity scores based on pubchem fingerprints

    Parameters:
    -----------
    basis_set (dict): dictionary of molecules with desired property in format {key:CID}. Molecules must be represented as pubchem CID or SMILES

    """ 
    assert isinstance(basis_set, dict), 'basis set must be a dictionary with Pubchem CIDs as values'
    assert utils.is_integery(list(basis_set.values())[0]), 'Basis set values must be Pubchem CIDs' 

    url_dict = {}
    for key in list(basis_set.keys()):
        cid = basis_set[key]
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/{cid}/cids/JSON?Threshold={threshold}&MaxRecords={max_records}'
        url_dict[key] = url
    
    retriever = qpc.pubchemQuery()
    similarity_responses = retriever.run_queries(url_dict)

    similarities = {}
    for key, value in similarity_responses.items():
        if value == 'FAILED':
            similarities[key] = value
        else:
            similarities[key] = value.json()['IdentifierList']['CID']

    return similarities

def get_molecule_properties(molecules: dict, properties: list):
    """
    Get the desired properties from pubchem for the molecules in molecules dictionary

    Parameters:
    -----------
    molecules: dictionary, values are pubchem CIDs
    properties: list of strings, properties to get from pubchem. Match pubchem property names, can be found here: 
    https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest#section=Compound-Property-Tables

    """
    assert isinstance(molecules, dict), 'basis set must be a dictionary with Pubchem CIDs as values'
    assert utils.is_integery(list(molecules.values())[0]), 'Basis set values must be Pubchem CIDs' 

    url_dict = {}
    for key in list(molecules.keys()):
        cid = molecules[key]
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/{",".join(properties)}/JSON'
        url_dict[key] = url
    
    retriever = qpc.pubchemQuery()
    property_responses = retriever.run_queries(url_dict)

    property_dict = {}

    for key, value in property_responses.items():
        if value == 'FAILED':
            property_dict[key] = value
        else:
            property_dict[key] = value.json()['PropertyTable']['Properties']

    return property_dict
