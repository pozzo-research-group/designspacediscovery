# functions related to running 2d similarity search
import designspacediscovery.querypubchem as qpc
import designspacediscovery.utils as utils
import requests
import pickle
import warnings

def find_similar_molecules(basis_set: dict,
                           threshold=90,
                           max_records=5000,
                           representation='CID') -> dict:
    """
    Find similar molecules using the pubchem fast2dsimilarity api

    Uses tanimoto similarity scores based on pubchem fingerprints

    Parameters:
    -----------
    basis_set (dict): dictionary of molecules with desired property in format {key:CID}. Molecules must be represented as pubchem CID or SMILES

    """
    assert isinstance(
        basis_set,
        dict), 'basis set must be a dictionary with Pubchem CIDs as values'
    assert utils.is_integery(list(
        basis_set.values())[0]), 'Basis set values must be Pubchem CIDs'

    url_dict = {}
    for key in list(basis_set.keys()):
        cid = basis_set[key]
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/{cid}/cids/JSON?Threshold={threshold}&MaxRecords={max_records}'
        url_dict[key] = url

    retriever = qpc.pubchemQuery()
    similarity_responses = retriever.run_queries(url_dict)

    similarities = {}
    for key, value in similarity_responses.items():
        #print(key, value)
        if value == 'FAILED':
            similarities[key] = value
        else:
            try:
                similarities[key] = value.json()['IdentifierList']['CID']
            except:
                print('exception on', key, value)
                similarities[key] = 'FAILED'

    return similarities


def get_molecule_properties(molecules: list, properties: list):
    """
    Get the desired properties from pubchem for the molecules in molecules dictionary

    Parameters:
    -----------
    molecules: list of pubchem CIDs
    properties: list of strings, properties to get from pubchem. Match pubchem property names, can be found here: 
    https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest#section=Compound-Property-Tables

    Returns:
    --------
    properties: dict with structure {key:{'CID':cid, 'Property1:prop1value, ...}}

    """
    assert isinstance(
        molecules, list), 'basis set must be a list of pubchem cids'
    assert utils.is_integery(molecules[0]), 'Basis set values must be Pubchem CIDs'

    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/{",".join(properties)}/JSON'
  
    retriever = qpc.pubchemQuery()
    property_responses = retriever.batch_queries(molecules, url)

    properties_list = []

    for resp in property_responses:
        if resp == "FAILED":
            pass
        elif not isinstance(resp, requests.models.Response):
            warnings.warn('Encountered unexpected response value in properties responses: ', value)
        elif resp.status_code !=200:
            warnings.warn('Bad response code for response: ', resp.status_code)
        else: 
            try:
                properties_list.extend(resp.json()['PropertyTable']['Properties'])
            except:
                warnings.warn('Error parsing json from otherwise valid response')
                print(resp.content)

    return properties_list

def get_pubchem_vendor_status(molecules, cache_params = {'cache':True, 'cache_fp':'.', 'cache_name':'vendor_cache'}):
    """
    Determine if a molecule is purchaseable based on pubchem vendors being present

    Uses serial pubchem queries. Expect to wait. 
    """
    url_dict = {}
    for key, cid in molecules.items():
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/{cid}/JSON'
        url_dict[key] = url

    retriever = qpc.pubchemQuery()
    vendor_responses = retriever.run_queries(url_dict, cache = cache_params['cache'], cache_fp = cache_params['cache_fp'], cache_name = cache_params['cache_name'])


    # cache the final result as a dict, this is a delicate fragile payload:
    if cache_params['cache']:
        with open(f'{cache_params["cache_fp"]}/{cache_params["cache_name"]}_final_vendors_responses.pkl', 'wb') as f:
            pickle.dump(vendor_responses, f)
    
    vendor_status = {}

    for key, value in vendor_responses.items():
        if value == "FAILED":
            vendor_status[key] = False
        elif not isinstance(value, requests.models.Response):
            vendor_status[key] = False
        elif value.status_code !=200:
            vendor_status[key] = False
        else:
            vendors = None
            # the response json is structured so that ['source cats']['cats'] is a list of dictionaries, one for each contributor category. Find the chemical vendor one, if it exists, and check how long it is  
            try:
                for cat in value.json()['SourceCategories']['Categories']:
                    if cat['Category'] == "Chemical Vendors":
                        vendors = cat
                        break
            except:
                pass

            if vendors is not None:
                if len(vendors) > 0:
                    vendor_status[key] = True
                else:
                    vendor_status[key] = False
            else: 
                vendor_status[key] = False
                
    return vendor_status
                            