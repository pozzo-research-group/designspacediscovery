import pytest
from designspacediscovery.vendors import is_purchaseable
from unittest.mock import patch


@patch('designspacediscovery.similaritysearch.get_pubchem_vendor_status')
def test_is_purchaseable_pubchem(mock_vendor_status):

    test_basis = {'1': 445639, '2': 18397614}
    test_pubchem_result = {'1': True, '2': False}

    mock_vendor_status.return_value = test_pubchem_result

    purchaseable = is_purchaseable(test_basis, criteria='pubchem')

    assert list(purchaseable.values()) == [True, False]

    return None
    # 1. test pubchem results


#@patch('designspacediscovery.similaritysearch.get_pubchem_vendor_status')
def test_is_purchaseable_bloom():

    test_basis = {'1': 445639, '2': 18397614}

    purchaseable = is_purchaseable(test_basis, criteria='bloom')

    assert list(purchaseable.values()) == [False, False]

    return None


    # 1. test pubchem results
@patch('designspacediscovery.similaritysearch.get_pubchem_vendor_status')
def test_is_purchaseable_either(mock_vendor_status):

    test_basis = {'1': 445639, '2': 18397614}
    test_pubchem_result = {'1': True, '2': False}

    mock_vendor_status.return_value = test_pubchem_result

    purchaseable = is_purchaseable(test_basis, criteria='either')

    assert list(purchaseable.values()) == [True, False]

    return None


@patch('designspacediscovery.similaritysearch.get_pubchem_vendor_status')
def test_is_purchaseable_both(mock_vendor_status):

    test_basis = {'1': 445639, '2': 18397614}
    test_pubchem_result = {'1': True, '2': False}

    mock_vendor_status.return_value = test_pubchem_result

    purchaseable = is_purchaseable(test_basis, criteria='both')

    assert list(purchaseable.values()) == [False, False]

    return None

    # 2. test bloom process
    # 3. Test both process
    # 4. test either process
