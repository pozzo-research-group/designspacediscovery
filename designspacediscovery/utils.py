import itertools
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
import tqdm
import numpy as np
import re
import contextlib
import sys


def is_integery(val):
    """
    Figure out if something could be a pubchem id
    """
    if (not isinstance(val, str)) and (not (isinstance(val, int) or isinstance(val, np.int64))):
        return False
    else:
        if isinstance(val, str):
            try:
                int(val)
                return True
            except:
                return False
        elif isinstance(val, int) or isinstance(val, np.int64):
            return True
        else:
            return False


def chunked_iterable(iterable, size):
    """
    Yields a chunk of size 'size' from the iterable. 
    """
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk


def fp_list_from_smiles_list(smiles_list, n_bits=2048):
    """
    Takes a list of smiles, returns a list of morgan fingerprints (rdkit). np.nan in list if failure.
    """
    fp_list = []
    for smiles in tqdm.tqdm(smiles_list, file=sys.stdout):
        with nostdout():
            #print('iterating')
            mol = Chem.MolFromSmiles(smiles)
            if mol == None:  #added this in to skip None as they returned sometimes in the line before
                fp_list.append(np.nan)
            else:
                fp_list.append(fp_as_array(mol, n_bits))
    return fp_list


def fp_as_array(mol, n_bits=2048):
    """
    get morgan fingerprint as array with rdkit
    """
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
    arr = np.zeros((1, ), int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def load_pubchem_conversion(fp):
    """
    Load the tab-separated id -> cid text file from pubchem identifier exchange
    """
    cid_dict = {}
    with open(fp, 'r') as f:
        for line in f:
            tci, cid = line.split('\t')
            cid = re.sub('[^0-9]', '', cid)
            if len(cid) > 0:
                cid_dict[tci] = cid
    return cid_dict


### context manager to stop rdkit writes from messing up progress bars


@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    save_stderr = sys.stderr
    sys.stdout = DummyFile(sys.stdout)
    sys.stderr = DummyFile(sys.stdout)
    try:
        yield
    finally:
        sys.stdout = save_stdout
        sys.stderr = save_stderr


class DummyFile():
    file = None

    def __init__(self, file):
        self.file = file

    def write(self, x):
        if len(x.rstrip()) > 0:
            tqdm.tqdm.write(x, file=self.file)
