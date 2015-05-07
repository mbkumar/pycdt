from pymatgen.matproj.rest import MPRester
from pymatgen.electronic_structure.bandstructure import BandStructure

def get_vbm(mpid, mapi_key=None):

    """
    Returns the valence band maxiumum (float) of the structure with
    MP-ID mpid.

    Args:
        mpid (str): MP-ID for which the valence band maximum is to
            be fetched from the Materials Project database
        mapi_key: Materials API key to access database
    """

    if not mapi_key:
        with MPRester() as mp:
            bs = mp.get_bandstructure_by_material_id(mpid)
    else:
        with MPRester(mapi_key) as mp:
            bs = mp.get_bandstructure_by_material_id(mpid)

    if bs is None:
        raise ValueError("Could not fetch band structure!")

    vbm_dict = bs.get_vbm()
    if vbm_dict is None:
        raise ValueError("No entry for valence band maximum!")

    return vbm_dict["energy"]


