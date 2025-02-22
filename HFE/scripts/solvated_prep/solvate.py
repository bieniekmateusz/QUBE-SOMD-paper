
# coding: utf-8

# Author: Julien Michel<br>
# Email:&nbsp;&nbsp; julien.michel@ed.ac.uk
# 
# # Solvate
# 
# Based on Molecular Setup from Lester Hedges. Solvates an input with a chosen water model.

# In[1]:


import BioSimSpace as BSS


# In[2]:


node = BSS.Gateway.Node("A node to solvate a molecule ready for molecular simulation with AMBER and GROMACS.")


# In[3]:


node.addAuthor(name="Julien Michel", email="julien.michel@ed.ac.uk", affiliation="University of Edinburgh")
node.setLicense("GPLv3")


# In[4]:
node.addInput("box_dim", BSS.Gateway.Float(help="The dimension of the one side of the box."))


node.addInput("input", BSS.Gateway.FileSet(help="A topology and coordinates file"))

node.addInput("water", BSS.Gateway.String(help="The name of the water model to use for solvation.",
                                          allowed=BSS.Solvent.waterModels(), default="tip3p"))

node.addInput("ion_conc", BSS.Gateway.Float(help="The ionic concentration in mol/litre.",
                                            minimum=0, maximum=1, default=0))

node.addInput("output", BSS.Gateway.String(help="The root name of the solvated output files."))


# In[5]:


node.addOutput("nodeoutput", BSS.Gateway.FileSet(help="The parameterised and solvated molecular system in AMBER and GROMACS format."))


# In[6]:


node.showControls()


# In[7]:


system = BSS.IO.readMolecules(node.getInput("input"))


# In[8]:

box_solv = [node.getInput("box_dim")*BSS.Units.Length.angstrom, node.getInput("box_dim")*BSS.Units.Length.angstrom, node.getInput("box_dim")*BSS.Units.Length.angstrom]
system = BSS.Solvent.solvate(node.getInput("water"), molecule=system,
                                                     box = box_solv,
                                                     ion_conc=node.getInput("ion_conc"))


# In[9]:


node.setOutput("nodeoutput", BSS.IO.saveMolecules(node.getInput("output"), system, ["pdb", "prm7", "rst7", "gro87", "grotop"]))


# In[10]:


node.validate()


# In[ ]:


#BSS.IO.saveMolecules(node.getInput("output"), system, ["prm7","rst7"])
