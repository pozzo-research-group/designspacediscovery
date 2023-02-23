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
@patch.object(designspacediscovery.querypubchem.pubchemQuery, 'run_queries')
@patch.object(requests.Response, 'json')
def test_get_molecule_properties_validurlcall_parseresponse(
        mock_response, mock_run_queries):
    # make sure a mock response is processed correctly
    correct_url_0 = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/445639/property/MolecularFormula,MolecularWeight/JSON'
    test_basis = {'1': 445639}
    test_properties = ['MolecularFormula', 'MolecularWeight']
    mock_response.json.return_value = {'IdentifierList': {'CID': [445639]}}
    mock_run_queries.return_value = {'1': mock_response}

    results = sim.find_similar_molecules(test_basis,
                                         max_records=1,
                                         threshold=99)

    assert list(mock_run_queries.call_args.args[0].values()
                )[0] == correct_url_0, 'URL passed to pubchem query malformed'

    assert results['1'] == [445639], 'Error in parsing of pubchem reponses'
