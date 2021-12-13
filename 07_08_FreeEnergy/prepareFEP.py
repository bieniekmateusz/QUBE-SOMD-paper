#!/usr/bin/env python
# coding: utf-8

# Author: Julien Michel
# 
# email: julien.michel@ed.ac.uk

# # PrepareFEP
# Loads a pair of input files, perform mapping between the first molecule of each input. Write down input files for a SOMD FEP calculation.

import BioSimSpace as BSS
import os
from Sire.Mol import AtomIdx
import shutil



def writeLog(ligA, ligB, mapping, output):
    """ Human readable report on atoms used for the mapping."""
    atoms_in_A = list(mapping.keys())
    stream = open(output,'w')
    atAdone = []
    atBdone= []
    for atAidx in atoms_in_A:
        # breakpoint()
        # import ipdb ; ipdb.set_trace()
        atA = ligA._getSireObject().select(AtomIdx(atAidx))
        atB = ligB._getSireObject().select(AtomIdx(mapping[atAidx]))
        stream.write("%s %s --> %s %s\n" % (atA.index(), atA.name(),atB.index(), atB.name()))
        atAdone.append(atA)
        atBdone.append(atB)
    for atom in ligA._getSireObject().atoms():
        if atom in atAdone:
            continue
        stream.write("%s %s --> dummy\n" % (atom.index(), atom.name()))
    for atom in ligB._getSireObject().atoms():
        if atom in atBdone:
            continue
        stream.write("dummy --> %s %s\n" % (atom.index(), atom.name()))
    stream.close()


def loadMapping(mapping_file):
    """Parse a text file that specifies mappings between atomic indices in input1 --> atoms in input2"""
    stream = open(mapping_file,'r')
    buffer = stream.readlines()
    stream.close()
    mapping = {}
    for line in buffer:
        if line.startswith("#"):
            continue
        elems = line.split(",")
        idx1 = int(elems[0])
        idx2 = int(elems[1])
        mapping[ AtomIdx(idx1)] = AtomIdx(idx2)
    
    return mapping


node = BSS.Gateway.Node("A node to generate input files for a SOMD relative free energy calculation.")

node.addAuthor(name="Julien Michel", email="julien.michel@ed.ac.uk", affiliation="University of Edinburgh")
node.setLicense("GPLv3")

node.addInput("input1", BSS.Gateway.FileSet(help="A topology and coordinates file"))
node.addInput("input2", BSS.Gateway.FileSet(help="A topology and coordinates file"))
node.addInput("prematch", BSS.Gateway.String(help="list of atom indices that are matched between input2 and input1. Syntax is of the format 1-3,4-8,9-11... Ignored if a mapping is provided", default=""))
node.addInput("mapping", BSS.Gateway.File(help="csv file that contains atom indices in input1 mapped ot atom indices in input2", optional=True))
node.addInput("output", BSS.Gateway.String(help="The root name for the files describing the perturbation input1->input2."))

node.addOutput("nodeoutput", BSS.Gateway.FileSet(help="SOMD input files for a perturbation of input1->input2."))

node.showControls()

do_mapping = True
custom_mapping = node.getInput("mapping")
if custom_mapping is not None:
    do_mapping = False
    mapping = loadMapping(custom_mapping)

# Optional input, dictionary of Atom indices that should be matched in the search. 
prematch = {}
prematchstring = node.getInput("prematch")
if len(prematchstring) > 0: 
    entries = prematchstring.split(",")
    for entry in entries:
        idxA, idxB = entry.split("-")
        prematch[ AtomIdx( int(idxA)) ] = AtomIdx( int(idxB) )

# Load system 1
system1 = BSS.IO.readMolecules(node.getInput("input1"))
# Load system 2
system2 = BSS.IO.readMolecules(node.getInput("input2"))

# We assume the molecules to perturb are the first molecules in each system
lig1 = system1.getMolecules()[0]
lig2 = system2.getMolecules()[0]

if do_mapping:
    # Return a maximum of 10 matches, scored by RMSD and sorted from best to worst.
    mappings, scores = BSS.Align.matchAtoms(lig1, lig2, matches=10, prematch=prematch, return_scores=True, scoring_function="RMSDalign", timeout=10*BSS.Units.Time.second)
    # We retain the top mapping
    mapping = mappings[0]

#print (mapping)
#for x in range(0,len(mappings)):
#    print (mappings[x], scores[x])

inverted_mapping = dict([[v,k] for k,v in mapping.items()])

# Align lig2 to lig1 based on the best mapping (inverted). The molecule is aligned based
# on a root mean squared displacement fit to find the optimal translation vector
# (as opposed to merely taking the difference of centroids).
lig2 = BSS.Align.rmsdAlign(lig2, lig1, inverted_mapping)
# Merge the two ligands based on the mapping.
merged = BSS.Align.merge(lig1, lig2, mapping)
# Create a composite system
system1.removeMolecules(lig1)
system1.addMolecules(merged)

root = node.getInput("output")
mapping_str = "%s.mapping" % root

# Log the mapping used
writeLog(lig1, lig2, mapping, output=mapping_str)
BSS.IO.saveMolecules(f"{root}.mergeat0.pdb", merged, "pdb", { "coordinates" : "coordinates0" })
# BSS.IO.saveMolecules("merged_at_lam0.pdb", merged, "pdb", { "coordinates" : "coordinates0" , "element": "element0" })
# Generate package specific input
protocol = BSS.Protocol.FreeEnergy(runtime = 2*BSS.Units.Time.nanosecond, num_lam=3)
#protocol = BSS.Protocol.FreeEnergy(runtime = 2*BSS.Units.Time.femtosecond, num_lam=3)
process = BSS.Process.Somd(system1, protocol)
process.getOutput() # creates the .zip? 
shutil.move('somd_output.zip', f'{root}.zip')

# we copy them directl from _work_dir ?
# cmd = "unzip -o somd_output.zip"
# os.system(cmd)

root = node.getInput("output")
pert = "%s.pert" % root
prm7 = "%s.prm7" % root
rst7 = "%s.rst7" % root

# fixme - shutil? 
cmd = f"""  mv {process._work_dir}/somd.pert {pert} ; 
            mv {process._work_dir}/somd.prm7 {prm7} ; 
            mv {process._work_dir}/somd.rst7 {rst7} ; 
       """
os.system(cmd)


node.setOutput("nodeoutput",[pert, prm7, rst7, mapping_str])
node.validate()