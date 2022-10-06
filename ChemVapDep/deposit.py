#!/usr/bin/env python

import os
import sys
import numpy as np
import shutil as sh
import gromacs as gmx
import argparse as arg
import MDAnalysis as mda

import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(
                formatter_class=arg.ArgumentDefaultsHelpFormatter)

    #
    # Input Options
    #
    inp = parser.add_argument_group("Input Data")

    inp.add_argument(
            '-s',
            '--surf',
            type=str,
            required=True,
            dest='SurfFile',
            help='''Surface coordinates.'''
        )

    inp.add_argument(
            '-m',
            '--mol',
            type=str,
            required=True,
            dest='MolFile',
            help='''Molecule coordinates.'''
        )

    inp.add_argument(
            '-a',
            '--ax',
            type=str,
            dest='DepAx',
            required=True,
            choices=[ "x", "y", "z" ],
            help='''Deposition axis.'''
        )

    inp.add_argument(
            '-mdp',
            '--mdp',
            type=str,
            required=True,
            dest='MDPFile',
            help='''Gromacs MDP template.'''
        )

    inp.add_argument(
            '-p',
            '--top',
            type=str,
            required=True,
            dest='TopFile',
            help='''Initial topology.'''
        )

    inp.add_argument(
            '-sp',
            '--surfitp',
            type=str,
            required=True,
            dest='SurfItpFile',
            help='''Surface Include Topology File.'''
        )

    inp.add_argument(
            '-sr',
            '--surfres',
            type=str,
            required=True,
            dest='SurfResName',
            help='''Surface Residue Name.'''
        )

    inp.add_argument(
            '-mp',
            '--molitp',
            type=str,
            required=True,
            dest='MolItpFile',
            help='''Molecule Include Topology File.'''
        )

    inp.add_argument(
            '-mr',
            '--molres',
            type=str,
            required=True,
            dest='MolResName',
            help='''Molecule Residue Name.'''
        )

    inp.add_argument(
            '-mt',
            '--maxtries',
            type=int,
            default=10000,
            dest='MaxTries',
            help='''Maximum number of attempts to add a new molecule.'''
        )

    inp.add_argument(
            '-n',
            '--nmol',
            type=int,
            default=216,
            dest='NMol',
            help='''Number of molecules to deposit.'''
        )

    args = parser.parse_args()
    Opts = vars(args)

    return Opts


def deposit_molecules(**Opts):

    cwd = os.getcwd()
    counter = 1

    while counter <= Opts['NMol']:

        # print()
        # print("----------------")
        # print("Starting loop %d" % counter)
        # print("----------------")
        # print()
        # Create folder for current molecule
        currentdir = "Mol_%03d" % counter
        currentdir = os.path.join(cwd, currentdir)
        os.makedirs(currentdir)

        # Different source of files if first loop or simulation in course
        if counter == 1:
            old = counter - 1
            olddir = cwd
        else:
            old = counter - 1
            olddir = os.path.join(cwd, "Mol_%03d" % old)

        # Copy necessary file to the folder
        files = [
            Opts['SurfFile'],
            Opts['MolFile'],
            Opts['TopFile'],
            Opts['MDPFile'],
            Opts['SurfItpFile'],
            Opts['MolItpFile']
        ]
        files = [ os.path.join(olddir, f) for f in files ]
        [ sh.copy(f, currentdir) for f in files ]

        # Move to the current directory to run gromacs
        os.chdir(currentdir)

        # Rename previous surface coordinates
        oldsurf = os.path.join(currentdir, Opts['SurfFile'])
        basename = os.path.join(currentdir, os.path.splitext(Opts['SurfFile'])[0])
        oldsurf = sh.move(oldsurf, basename + "_%03d.gro" % old)
        newsurf = basename + "_%03d_ini.gro" % counter

        # Rename previous topology
        oldtop = os.path.join(currentdir, Opts['TopFile'])
        oldtop = sh.copy(oldtop, oldtop + ".old")

        # Select atoms with the highest coordinate along the deposition axis
        # Conversion to nm, as it is the gromacs input
        if Opts["DepAx"] == "x":
            idx = 0
        elif Opts["DepAx"] == "y":
            idx = 1
        else:
            idx = 2

        molres = Opts["MolResName"]
        surfres = Opts["SurfResName"]
        u = mda.Universe(os.path.basename(oldsurf))
        box = u.dimensions

        # Error handling for first loop when no molecules are deposited
        # requesting the max from an empty array gives a ValueError
        try:
            mollow = u.select_atoms(f"resname {molres}").positions[:,idx].max() / 10.0
        except ValueError:
            mollow = -np.inf

        surflow = u.select_atoms(f"resname {surfres}").positions[:,idx].max() / 10.0

        # Define the region of space along the deposition axis where the new
        # molecule can be selected at the beginning of the simulation.
        # Let's say 10-15 A from the highest coordinate along the deposition
        # axis, but this is totally arbitrary and depends on the LJ curve.
        low = np.max([ mollow, surflow ]) + 1.0
        high = low + 0.5

        # Add one molecule to the simulation
        # Maybe can be improved to add more than one, but this can be easily
        # solved by providing an input file for the molecule to be deposited
        # that contains more than one residue.
        result = False
        attempts = 0

        # Here we try adding the molecule until we are successful.
        # This could be done with the -try option of gmx.insert_molecules, but
        # in anycase we would need error handling in case it fails after the
        # maximum number of tries.
        while not result:
            try:
                tlow = np.zeros(3)
                tlow[idx] = low
                thigh = box[:3]
                thigh[idx] = high
                pos = np.random.uniform(
                        size=3,
                        low=tlow,
                        high=thigh
                        # low=[ low + 1.0, 0.0, 0.0 ],
                        # high=[ high, 6.0, 6.0 ]
                    )
                with open("positions.dat", "w") as ip:
                    ip.write("%4.2f %4.2f %4.2f\n" % tuple(pos.tolist()))

                gmx.insert_molecules(
                        f=os.path.basename(oldsurf),
                        ci=Opts['MolFile'],
                        nmol=1,
                        rot="xyz",
                        ip="positions.dat",
                        o=os.path.basename(newsurf)
                    )

                # Make index file
                gmx.make_ndx(
                        input=('q'),
                        f=os.path.basename(newsurf),
                    )

                # Select atoms with low z coordinate to keep the surface in place
                ax = Opts["DepAx"]
                g = u.select_atoms(f"resname PCB and prop {ax} <= 4.0")
                with mda.selections.gromacs.SelectionWriter('index.ndx', mode='a') as ndx:
                    ndx.write(g, name='base')

                # Update topology adding one molecule
                topold = os.path.join(currentdir, Opts['TopFile'] + ".old")
                with open(topold) as top:
                    lines = top.readlines()

                res, prev = lines[-1].strip().split()[:2]
                lines[-1] = "%s                    %d" % ( res, int(prev) + 1 )

                topname = os.path.join(currentdir, Opts['TopFile'])
                with open(topname, "w") as top:
                    top.writelines(lines)

                # Grompp
                mdpfile = os.path.join(currentdir, Opts['MDPFile'])
                tprfile = "deposit_%03d.tpr" % counter
                gmx.grompp(
                        f=os.path.basename(mdpfile),
                        c=os.path.basename(newsurf),
                        p=os.path.basename(topname),
                        n="index.ndx",
                        o=tprfile,
                        maxwarn=1
                    )

                result = True
                attempts += 1

            # Every N failed attempts, the region of space to insert the
            # new molecule gets increased by 10 A. The program stops if we
            # are not able to add the molecule after a fixed number of tries.
            except:
                attempts += 1
                if attempts % (Opts["MaxTries"] // 50) == 0:
                    high += 1.0
                elif attempts > Opts["MaxTries"]:
                    break

        # Run MD
        tprfile_base = "deposit_%03d" % counter
        gmx.mdrun(
                deffnm=tprfile_base,
                v=True,
                nt=16
            )

        # Rename output coordinates for next step
        sh.copy("%s.gro" % tprfile_base, Opts['SurfFile'])

        # Update counter
        counter += 1

    return


def main():

    Opts = options()
    deposit_molecules(**Opts)

    return


if __name__ == '__main__':
    main()
