"""
Get the possible vendors for a chemical

2 methods: molbloom/zinc or pubchem vendor list

"""
import molbloom
import designspacediscovery.similaritysearch as sim
import designspacediscovery.utils as utils


def is_purchaseable(molecules, criteria, smiles_dict=None):
    """
    Lookup for purchase-ability of compounds
    
    Parameters:
    ----------
    molecules - dictionary with cids as values
    criteria (bloom, pubchem, either, both) - which source or combo of sources to use
    smiles_dict: optional, dictionary of smiles strings for the keys in molecules. Smiles strings required for bloom lookup. If bloom required and smiles_dict missing, will use pubchem to get smiles strings
    """

    assert isinstance(
        molecules,
        dict), "This function takes a dictionary of smiles strings as input"
    assert utils.is_integery(list(
        molecules.values())[0]), 'Dictionary keys should be smiles strings'
    assert criteria in [
        'pubchem', 'bloom', 'either', 'both'
    ], "Please select a criteria from 'pubchem', 'bloom', 'either', 'both'"

    if criteria in ['pubchem', 'either', 'both']:
        # run pubchem search
        pubchem_results = sim.get_pubchem_vendor_status(molecules)

    if criteria in ['bloom', 'either', 'both']:
        # if user did not supply a dictionary of smiles strings, go get it for them
        if smiles_dict is None:
            names = list(molecules.keys())
            smiles_list = sim.get_molecule_properties(
                molecules, properties=['CanonicalSMILES'])
            assert len(smiles_list) == len(
                names
            ), 'List of returned smiles is different than number of molecules'
        else:
            pass
        bloom_results = {
            name: molbloom.buy(v['CanonicalSMILES'])
            for name, v in zip(names, smiles_list)
        }

    if criteria == 'bloom':
        purchaseable = bloom_results
    elif criteria == 'pubchem':
        purchaseable = pubchem_results
    else:
        purchaseable = {}
        for key, cid in molecules.items():
            if criteria == 'both':
                purchaseable[key] = bloom_results[key] and pubchem_results[key]
            elif criteria == 'either':
                purchaseable[key] = bloom_results[key] or pubchem_results[key]

    return purchaseable
