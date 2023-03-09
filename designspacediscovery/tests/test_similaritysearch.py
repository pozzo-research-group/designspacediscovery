import requests
from unittest.mock import patch
import time
import designspacediscovery.similaritysearch as sim
import designspacediscovery.querypubchem
import pytest


def test_find_similar_molecules_inputchecking():
    good_dict = {'1': 445639, '2': 11005, '3': 5281, '4': 3893}
    bad_dict = {'1': 'aabbd112'}

    # make sure input checking works
    with pytest.raises(AssertionError):
        sim.find_similar_molecules(bad_dict)
    with pytest.raises(AssertionError):
        sim.find_similar_molecules([1234, 5332])


@patch.object(designspacediscovery.querypubchem.pubchemQuery, 'run_queries')
@patch.object(requests.Response, 'json')
def test_find_similar_molecules_validurlcall_parseresponse(
        mock_response, mock_run_queries):
    # make sure a mock response is processed correctly
    correct_url_0 = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/445639/cids/JSON?Threshold=99&MaxRecords=1'
    test_basis = {'1': 445639}
    mock_response.json.return_value = {'IdentifierList': {'CID': [445639]}}
    mock_run_queries.return_value = {'1': mock_response}

    results = sim.find_similar_molecules(test_basis,
                                         max_records=1,
                                         threshold=99)

    assert list(mock_run_queries.call_args.args[0].values()
                )[0] == correct_url_0, 'URL passed to pubchem query malformed'

    assert results['1'] == [445639], 'Error in parsing of pubchem reponses'


# make sure URL is valid
@patch.object(designspacediscovery.querypubchem.pubchemQuery, 'batch_queries')
@patch.object(requests.models.Response, 'json')
def test_get_molecule_properties_validurlcall_parseresponse(
        mock_response, mock_run_queries):
    # make sure a mock response is processed correctly
    correct_url_0 = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/MolecularFormula,MolecularWeight/JSON'
    test_basis = {'1': 445639}
    test_properties = ['MolecularFormula', 'MolecularWeight']
    mock_response.json.return_value = {
        'PropertyTable': {
            'Properties': [{
                'CID': 445639,
                'MolecularFormula': 'C18H34O2',
                'MolecularWeight': '282.5'
            }]
        }
    }
    mock_response.status_code = 200
    mock_run_queries.return_value = [mock_response]

    results = sim.get_molecule_properties(test_basis, test_properties)

    assert mock_run_queries.call_args.args[0][0] == test_basis[
        '1'], 'Wrong list passed to batch query'
    assert mock_run_queries.call_args.args[
        1] == correct_url_0, 'Wrong URL passed to batch query'

    assert list(
        results[0].keys())[0] == 'CID', 'Error in parsing of pubchem reponses'


@patch.object(designspacediscovery.querypubchem.pubchemQuery, 'run_queries')
@patch.object(requests.Response, 'json')
def test_get_pubchem_vendor_status(mock_response, mock_run_queries):
    test_basis = {'1': 445639}

    correct_url_0 = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/445639/JSON'
    mock_json = {
        'SourceCategories': {
            'RecordType':
            'CID',
            'RecordNumber':
            445639,
            'Categories': [{
                'Category':
                'Chemical Vendors',
                'URL':
                'https://www.ncbi.nlm.nih.gov/sites/entrez?cmd=search&db=pcsubstance&term=445639[StandardizedCID]+AND+%22Chemical+Vendors%22[SourceCategory]',
                'Sources': [{
                    'SID': 441561215,
                    'SourceName': 'Elsa Biotechnology',
                    'SourceURL': 'https://www.elsa-biotech.com',
                    'SourceDetail':
                    'https://pubchem.ncbi.nlm.nih.gov/source/Elsa%20Biotechnology',
                    'RegistryID': 'ELSA20-8625'
                }, {
                    'SID':
                    355047057,
                    'SourceName':
                    'Yuhao Chemical',
                    'SourceURL':
                    'http://www.chemyuhao.com',
                    'SourceDetail':
                    'https://pubchem.ncbi.nlm.nih.gov/source/Yuhao%20Chemical',
                    'RegistryID':
                    'HL1416',
                    'SourceRecordURL':
                    'http://www.chemyuhao.com/112-80-1.html'
                }]
            }]
        }
    }
    mock_response.json.return_value = mock_json
    mock_response.status_code = 200
    mock_run_queries.return_value = {'1': mock_response}

    # caching logic is not checked here due to issue with picking mock objects
    results = sim.get_pubchem_vendor_status(test_basis,
                                            cache_params={'cache': False})

    # test that retriever is called correctly
    assert list(mock_run_queries.call_args.args[0].values(
    ))[0] == correct_url_0, 'Incorrect url passed for pubchem vendor getting'

    # test that a mock result gets parsed correctly
    assert results['1'] == True, 'Incorrectly parsed pubchem vendor response'

    return None


@patch.object(designspacediscovery.querypubchem.pubchemQuery, 'run_queries')
@patch.object(requests.Response, 'json')
def test_get_pubchem_vendor_status_false(mock_response, mock_run_queries):
    test_basis = {'1': 445639}

    correct_url_0 = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/445639/JSON'
    mock_json = {
        'SourceCategories': {
            'RecordType':
            'CID',
            'RecordNumber':
            445639,
            'Categories': [{
                'Category':
                'Rainbows',
                'URL':
                'https://www.ncbi.nlm.nih.gov/sites/entrez?cmd=search&db=pcsubstance&term=445639[StandardizedCID]+AND+%22Chemical+Vendors%22[SourceCategory]',
                'Sources': [{
                    'SID': 441561215,
                    'SourceName': 'Elsa Biotechnology',
                    'SourceURL': 'https://www.elsa-biotech.com',
                    'SourceDetail':
                    'https://pubchem.ncbi.nlm.nih.gov/source/Elsa%20Biotechnology',
                    'RegistryID': 'ELSA20-8625'
                }, {
                    'SID':
                    355047057,
                    'SourceName':
                    'Yuhao Chemical',
                    'SourceURL':
                    'http://www.chemyuhao.com',
                    'SourceDetail':
                    'https://pubchem.ncbi.nlm.nih.gov/source/Yuhao%20Chemical',
                    'RegistryID':
                    'HL1416',
                    'SourceRecordURL':
                    'http://www.chemyuhao.com/112-80-1.html'
                }]
            }]
        }
    }
    mock_response.json.return_value = mock_json
    mock_response.status_code = 200
    mock_run_queries.return_value = {'1': mock_response}

    # caching logic is not checked here due to issue with picking mock objects
    results = sim.get_pubchem_vendor_status(test_basis,
                                            cache_params={'cache': False})

    # test that retriever is called correctly
    assert list(mock_run_queries.call_args.args[0].values(
    ))[0] == correct_url_0, 'Incorrect url passed for pubchem vendor getting'

    # test that a mock result gets parsed correctly
    assert results[
        '1'] == False, 'Incorrectly parsed pubchem vendor response for false example'

    return None


def test_substructure_search():
    acrylate_smarts = '[#6&D1]=[#6&D2]-[#6](=[#8])-[#8]'

    test_smiles = {
        'yes': 'C=CC[Si](CCCOC(=O)C=C)(CC=C)CC=C',
        'no': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O'
    }

    match_dict = sim.substructure_search(acrylate_smarts, test_smiles)

    assert match_dict['yes'] == True, 'match not found'
    assert match_dict['no'] == False, 'Not match flagged as match'
