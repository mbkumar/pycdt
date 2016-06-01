#!/usr/bin/env python

from future import divison

"""
Script file for end users to parse the calculations and get the 
defects data for further processing.
"""

from arg_parser import ArgParse
import argparser

from pymatgen.matproj.rest import MPRester
from pycdcd.utilities.parsing_calculations import  parse_defect_calculations, \
        get_vbm_bandgap, get_atomic_chempots, parse_dielectric_calculation
from pycdcd.corrections.defect_analyzer import  DefectAnalyzer
from pycdcd.utils.plotter import  DefectPlotter

def get_parsed_data(mpid, MAPI_key=None, calc_root_fldr=None):
    """
    Function to parse the defects data computed for mpid.
    If the root folder for calculation is not given, deduce
    it form mpid based on how the root folder is generated from mpid.
    """
    if not calc_root_fldr:
        calc_root_fldr = mpid
        
    parsed_data = parse_defect_calculations(calc_root_fldr)
    epsilon = parse_dielectric_calculation(calc_root_fldr)
    blk_entry = parsed_data['bulk_entry']
    parsed_defects = parsed_data['defects']
    for defect in parsed_defects:
        apply_correction(defect, blk_entry, epsilon) # Freysoldt correction
    vbm, bandgap = get_vbm_bandgap(mpid, MAPI_key)
    struct = MPRester().get_structure_by_material_id()
    mus = get_atomic_chempots(struct)
    da = DefectAnalyzer(blk_entry, vbm, mus, bandgap['energy'])
    for defect in parsed_defects:
        da.add_parsed_defect(defect)
    return da


def get_formation_energy_plot(mpid, MAPI_key=None, calc_root_fldr=None):
    da = get_parsed_data(mpid)
    plotter = DefectPlotter(da)
    form_en_plot = plotter.get_plot_form_energy()
    form_en_plot.savefig(mpid+'formation_energy.eps')

