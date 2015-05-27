#!/Users/nisse/code/pymatgen/python2.7_devSIAs/bin/python

from __future__ import division, print_function, unicode_literals

"""
A script with tools for computing formation energies
of charged point defects, supporting multiple correction
schemes.
"""

__author__ = "Nils E. R. Zimmermann, Bharat Medasani"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1"
__maintainer__ = "Nils E. R. Zimmermann"
__email__ = "nerz@lbl.gov"
__date__ = "May 1, 2015"

import argparse
#import os
#import json
#import glob
import math

import pymatgen
from pymatgen.matproj.rest import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.standard_transformations import SupercellTransformation
from pymatgen.analysis.bond_valence import BVAnalyzer
from pycdcd.core.defectsmaker import ChargedDefectsStructures
from pycdcd.utils.vasp import make_vasp_defect_files, \
        make_vasp_dielectric_files
from pycdcd.utils.parse_calculations import PostProcess


def print_error_message(err_str):
    print("\n=================================================================="
        "=============\n\nError: "+err_str)
    print("\n================================================================"
        "===============\n")


def generate_input_files(args):
    """
    Generates input files for VASP calculations that aim to determine
    formation energies of charged point defects by (possibly) applying
    correction terms (supported so-far: correction due to Freysoldt
    et al., Phys. Rev. Lett., 2009).
    The primitive unit cell is obtained from the MP ID provided during
    script call.  

    Args:
        args (Namespace): contains the parsed command-line arguments for
            this command.
    """

    # initialize variables
    mpid = args.mpid
    mapi_key = args.mapi_key
    nmax = args.nmax
    oxi_state = args.oxi_state
    oxi_range = args.oxi_range

    # error-checking
    if not mpid:
        print_error_message("No Materials Project ID (MP-ID) provided!")
        return
    if nmax <= 0:
        print_error_message("maximal number of atoms per supercell"
            " must be larger than zero!")
        return

    # get primitive unit cell
    if not mapi_key:
        with MPRester() as mp:
            prim_struct = mp.get_structure_by_material_id(mpid)
    else:
        with MPRester(mapi_key) as mp:
            prim_struct = mp.get_structure_by_material_id(mpid)

    # transform to conventional unit cell
    conv_struct = SpacegroupAnalyzer(
        prim_struct).get_conventional_standard_structure()

    make_vasp_dielectric_files(prim_struct)

    # manually set oxidation states if those were provided
    oxi_state_dict = {}
    if oxi_state:
        for i in range(len(oxi_state)):
            oxi_state_dict[oxi_state[i][0]] = int(oxi_state[i][1])
        if len(oxi_state_dict) != conv_struct.ntypesp:
            print_error_message("number of oxidation states"
                " provided does not match number of species in structure!")
            return

    # manually set oxidation-state ranges if those were provided
    oxi_range_dict = {}
    if oxi_range:
        for i in range(len(oxi_range)):
            oxi_range_dict[oxi_range[i][0]] = tuple([int(oxi_range[i][1]),
                int(oxi_range[i][2])])
        if len(oxi_range_dict) != conv_struct.ntypesp:
            print_error_message("number of distinct oxidation ranges"
                " provided does not match number of species in structure!")
            return

    # finally, generate VASP input files for defect calculations
    def_structs = ChargedDefectsStructures(conv_struct, 
            max_min_oxi=oxi_range_dict, oxi_states=oxi_state_dict, 
            cellmax=nmax)
    make_vasp_defect_files(def_structs.defects,
            conv_struct.composition.reduced_formula)


def parse_vasp_output(args):
    """
    Parses output files from VASP calculations that aim to determine
    formation energies of charged point defects by (possibly) applying
    correction terms (supported so-far: correction due to Freysoldt
    et al., Phys. Rev. Lett., 2009).

    Args:
        args (Namespace): contains the parsed command-line arguments for
            this command.
    """

    # initialize variables
    mpid = args.mpid
    mapi_key = args.mapi_key
    root_fldr = args.root_fldr

    # error-checking
    if not mpid:
        print_error_message("no Materials Project structure ID (MP-ID) provided!")
        return

    # get primitive unit cell
    if not mapi_key:
        with MPRester() as mp:
            prim_struct = mp.get_structure_by_material_id(mpid)
    else:
        with MPRester(mapi_key) as mp:
            prim_struct = mp.get_structure_by_material_id(mpid)

    output_dict = PostProcess(root_fldr, mpid, mapi_key).compile_all()
    #print(output_dict)
    #for pd in output_dict['defects']:
    #    print(pd.as_dict())

def main():
    parser = argparse.ArgumentParser(description="""
        pycdt is a script that generates vasp input files, parses vasp output
        files, and computes the formation energy of charged defects.
        This script works based on several sub-commands with their own options.
        To see the options for the sub-commands, type
        "pycdt sub-command -h".""",
        epilog="""
        Authors: N. E. R. Zimmermann, B. Medasani, D. Broberg, G. Hautier
        Version: {}
        Last updated: {}""".format(__version__, __date__))

    subparsers = parser.add_subparsers()
    mpid_string = "Materials Project id of the structure.\nFor more info on " \
        "Materials Project, please, visit www.materialsproject.org"
    mapi_string = "Your Materials Project REST API key.\nFor more info, " \
        "please, visit www.materialsproject.org/open"
    nmax_string = "Maximum number of atoms in supercell.\nThe default is" \
        "128.\nKeep in mind the number of atoms in the supercell may vary" \
        "from the provided number including the default."
    oxi_state_string = "Oxidation state for an element.\nTwo arguments" \
        " expected: the element type for which the oxidation state is" \
        " to be specified and the oxidation state (e.g., --oxi_state As -3)."
    oxi_range_string = "Oxidation range for an element.\nThree arguments" \
        " expected: the element type for which the oxidation state range is" \
        " to be specified as well as the lower and the upper limit of the" \
        " range (e.g., --oxi_range As -3 5)."
    root_fldr_string = "Path (relative or absolute) to directory" \
        " in which data of charged point-defect calculations for" \
        " a particular system are to be found.\n"

    parser_input_files = subparsers.add_parser("generate_input_files",
        help="Generates Vasp input files for charged point defects.")
    parser_input_files.add_argument("--mpid", type=str.lower, dest="mpid",
        help=mpid_string)
    parser_input_files.add_argument("--mapi_key", default=None,
        dest="mapi_key", help=mapi_string)
    parser_input_files.add_argument("--nmax", type=int, default=128,
        dest="nmax", help=nmax_string)
    parser_input_files.add_argument("--oxi_range", action='append', type=str,
        nargs=3, dest="oxi_range", help=oxi_range_string)
    parser_input_files.add_argument("--oxi_state", action='append', type=str,
        nargs=2, dest="oxi_state", help=oxi_state_string)
    parser_input_files.set_defaults(func=generate_input_files)

    parser_vasp_output = subparsers.add_parser("parse_vasp_output",
        help="Parses VASP output for calculation of formation energies of"
             " charged point defects.")
    parser_vasp_output.add_argument("--mpid", type=str.lower, dest="mpid",
        help=mpid_string)
    parser_vasp_output.add_argument("--mapi_key", default = None,
        dest="mapi_key", help=mapi_string)
    parser_vasp_output.add_argument("--dir", default = None,
        dest="root_fldr", help=root_fldr_string)
    parser_vasp_output.set_defaults(func=parse_vasp_output)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
