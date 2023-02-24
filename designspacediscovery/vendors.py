"""
Get the possible vendors for a chemical

2 methods: molbloom/zinc or pubchem vendor list

"""
import molbloom

def is_purchaseable(smiles_dict, criteria):
    """
    Lookup for purchasbility of compounds
    
    Parameters:
    ----------
    smiles (str) - dictionary with smiles strings as values
    criteria (bloom, pubchem, either, both) - which source or combo of sources to use
    """

    assert isinstance(smiles_dict, dict), "This function takes a dictionary of smiles strings as input"
    assert isinstance(list(smiles_dict.values())[0], str), 'Dictionary keys should be smiles strings'
    

    if criteria in ['pubchem', 'either', 'both']:
        # run pubchem search
        pubchem_results = pubchem_vendor_status(smiles_dict)
    
    if criteria in ['bloom', 'either', 'both']:
        bloom_results = {k:molbloom.buy(v) for k, v in smiles_dict.items()}

        