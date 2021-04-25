# This script can be used for any purpose without limitation subject  to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2020-01-21: created by S.B.Wiggin, the Cambridge Crystallographic Data Centre,
# and Aurelia Li, Adsorption and Advanced Materials Group (aam.ceb.cam.ac.uk),
# led by David Fairen-Jimenez from the Department of Chemical Engineering and
# Biotechnology, University of Cambridge.


"""
    Classify_MOFs_dimensionality.py  - performs two expansions of a polymeric network and produces minimum area bounding
    boxes. Comparison of the two bounding boxes gives the growth dimensions of the framework
"""

from ccdc.io import EntryReader
import numpy as np
import argparse
import csv
import os.path
import time


def generate_bounding_box(atoms):

    all_pts = np.array([[c for c in atom.coordinates] for atom in atoms])
    cov = np.cov(all_pts, rowvar=False)
    evals, evecs = np.linalg.eig(cov)

    lengths = np.sqrt(evals)
    lengths = lengths + 0.0001  # required when the box dimension is exceedingly small and tiny changes effect the ratio
    lengths.sort()

    return lengths


def dimensionality(entry):
    """
    Calculates the dimensionality of the crystal.
    Grows 2 instances of the polymeric unit (four cycles and seven cycles) and calculates the change in size
    The number of cycles needs to grow a representative part of the framework, but more cycles will take longer
    """

    t0 = time.time()
    shell1 = entry.crystal.polymer_expansion(repetitions=4)
    shell2 = entry.crystal.polymer_expansion(repetitions=7)

    t1 = time.time()
    time_taken = str(t1-t0)
    print('total time elapsed for polymer expansion is % s' % time_taken)

    lengths1 = generate_bounding_box(shell1.atoms)
    print('number of atoms in first expansion ' + str(len(shell1.atoms)))
    s1, m1, l1 = lengths1

    lengths2 = generate_bounding_box(shell2.atoms)
    print('number of atoms in second expansion ' + str(len(shell2.atoms)))
    s2, m2, l2 = lengths2

    ratios = [s2 / s1, m2 / m1, l2 / l1]
    print('ratio of box dimensions from first and second expansions: ' + str(ratios))
    thresh = 1.15
    ndims = [r > thresh for r in ratios].count(True)
    return ndims


def analyse_structures(user_gcd_input, user_csv_output):

    if len(os.path.splitext(user_csv_output)[1]) == 0:
        user_csv_output += ".csv"

    with open(user_csv_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(('Refcode', 'dimensionality', 'number in gcd file'))

        csd_reader = EntryReader(user_gcd_input, 'CSD')

        t2 = time.time()
        n_structures = 0
        n_mof = 0
        n_non_mof = 0

        for entry in csd_reader:
            print('CSD entry: ' + str(entry.identifier))
            n_structures += 1  # quick counter
            count_polymers = 0
            for component in entry.molecule.components:
                if component.is_polymeric:
                    count_polymers += 1
            if count_polymers > 1:
                print('multiple polymer units present')
            if entry.molecule.heaviest_component.is_polymeric:
                n_mof += 1
                framework = entry.molecule.heaviest_component

                framework.remove_hydrogens()  # next steps fail if any atoms in the unit do not have coordinates

                entry.crystal.molecule = framework

                fig = dimensionality(entry)

                if fig == 0:
                    dimension = '0D non-MOF'
                elif fig == 1:
                    dimension = '1D chain'
                elif fig == 2:
                    dimension = '2D sheet'
                elif fig == 3:
                    dimension = '3D framework'
            else:
                n_non_mof += 1
                dimension = 'no polymeric bonds detected'

            print('Framework dimensions for CSD entry % s: % s \n' % (entry.identifier, dimension))
            writer.writerow((entry.identifier, dimension, n_structures))
            f.flush()

        print('Total MOF subset size is: % d' % n_structures)
        print('Entries recognised as polyermic is: % d' % n_mof)
        print('Entries not recognised as polymeric (and ignored) is: % d' % n_non_mof)

        t3 = time.time()
        overall_time_taken = str(t3 - t2)
        print('total time elapsed for script % s' % overall_time_taken)
        f.close()


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', help='CSD refcode list (.gcd file) filename')
    parser.add_argument('-o', '--output', help='Results CSV filename')

    return parser.parse_args()


def main():
    args = get_args()

    if args.input is None:
        args.input = input("Enter filepath for CSD refcode list (.gcd file)" '\n')

    if args.output is None:
        args.output = input("Enter filename for results file" '\n')

    analyse_structures(args.input, args.output)


if __name__ == '__main__':
    main()
