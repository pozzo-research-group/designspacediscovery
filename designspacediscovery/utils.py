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
