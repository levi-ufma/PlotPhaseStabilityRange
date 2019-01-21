#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 16:53:08 2018

@author: jherfson
"""

"""
This is a basic example of how to create, plot, and analyze OPEN Phase Diagrams using the pymatgen
codebase and Materials Project database. To run this example, you should:
* have pymatgen (www.pymatgen.org) installed along with matplotlib
* obtain a Materials Project API key (https://www.materialsproject.org/open)
* paste that API key in the MAPI_KEY variable below, e.g. MAPI_KEY = "foobar1234"
For citation, see https://www.materialsproject.org/citing
For the accompanying comic book, see http://www.hackingmaterials.com/pdcomic
"""

# from __future__ import print_function
from pymatgen import MPRester, Element
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram, \
    PhaseDiagram, PDPlotter


#from graf.mu_to_T import print_T_corresponding_to_mu_equals
# import mu_to_T
from mu_to_temp import mu_to_temperature

imprimir = mu_to_temperature()


def plot_pd(pd, show_unstable=False):
    plotter = PDPlotter(pd, show_unstable=show_unstable)
    # plotter = PDPlotter(pd, show_unstable=0.03)

    plotter.show()
    # plotter.write_image("{}.png".format('-'.join(system)), "png")  # save figure


def analyze_pd(pd):
    print('Stable Entries (formula, materials_id)\n--------')
    for e in pd.stable_entries:
        print(e.composition.reduced_formula, e.entry_id)

    print('\nUnstable Entries (formula, materials_id, e_above_hull (eV/atom), decomposes_to)\n--------')
    for e in pd.unstable_entries:
        decomp, e_above_hull = pd.get_decomp_and_e_above_hull(e)
        pretty_decomp = [("{}:{}".format(k.composition.reduced_formula, k.entry_id), round(v, 2)) for k, v in
                         decomp.items()]
        print(e.composition.reduced_formula, e.entry_id, "%.3f" % e_above_hull, pretty_decomp)


def get_gcpd_data(gcpdobj=None):
    if gcpdobj is None:
        return []

    #  potential = [v for el, v in gcpdobj.chempots.items()]
    phases = ", ".join([entry.name for entry in gcpdobj.stable_entries])

    return [phases]


if __name__ == "__main__":
    MAPI_KEY = "Fvlb5EsNq71JxDy3"  # You must change this to your Materials API key! (or set MAPI_KEY env variable)
    system = ["Li", "B", "O"]  # system we want to get open PD for
    # system = ["Li", "Fe", "P", "O"]  # alternate system example

    open_elements_specific = None  # e.g., {Element("O"): 0} where 0 is the specific chemical potential
    open_element_all = Element("O")  # plot a series of open phase diagrams at critical chem pots with this element open

    mpr = MPRester(MAPI_KEY)  # object for connecting to MP Rest interface

    # get data
    entries = mpr.get_entries_in_chemsys(system, compatible_only=True)

    if open_elements_specific:
        gcpd = GrandPotentialPhaseDiagram(entries, open_elements_specific)
        plot_pd(gcpd, False)
        analyze_pd(gcpd)

    if open_element_all:
        pd = PhaseDiagram(entries)
        chempots = pd.get_transition_chempots(open_element_all)
        all_gcpds = list()
        toplot = []
        arquivo = open("dados.txt", "w")
        for idx in range(len(chempots)):
            if idx == len(chempots) - 1:
                avgchempot = chempots[idx] - 0.1
            else:
                avgchempot = 0.5 * (chempots[idx] + chempots[idx + 1])
            gcpd = GrandPotentialPhaseDiagram(entries, {open_element_all: avgchempot}, pd.elements)

            min_chempot = None if idx == len(chempots) - 1 else chempots[idx + 1]
            max_chempot = chempots[idx]
            #  print("Chempot range for diagram {} is: {} to {}".format(idx, min_chempot, max_chempot))

            toplot.append(get_gcpd_data(gcpd))

            print(get_gcpd_data(gcpd), min_chempot, max_chempot)
            #printf1 = imprimir.print_temperature_corresponding_to_mu_equals(max_chempot)
            #printf2 = imprimir.print_temperature_corresponding_to_mu_equals(min_chempot)
            #print(printf1, printf2)

            #print('===========================================')


            plot_pd(gcpd, False)









