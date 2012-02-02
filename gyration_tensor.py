#!/usr/bin/env python

# Copyright Notice
# ================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you are
# free to use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# -----------------------------------------------------------------------------------
# This PyMOL Plugin Gyration Tensor is
# Copyright (C) 2012 by Venkatramanan Krishnamani <venks@andrew.cmu.edu>
# 
#                        All Rights Reserved
# 
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
# 
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
#

###############################################################################
#                        Gyration Tensor 1.0                                  #
###############################################################################

import Tkinter
import tkFileDialog, tkMessageBox
import numpy
from pymol import cmd
from pymol.cgo import *

dialog = Tkinter.Tk()
dialog.withdraw()

def __init__(self):
        self.menuBar.addmenuitem('Plugin', 'command','Gyration Tensor',label = 'Gyration Tensor',command = lambda s=self : open_pdb(s))
        
def open_pdb(self):
        myFormatsPDB = [('Protein Data Bank','*.pdb'), ('MDL mol','*.mol'), ('PyMol Session File','*.pse')]
        try:
                self.pdbFile = tkFileDialog.askopenfile(parent=dialog,mode='rb',filetypes=myFormatsPDB, title='Choose PDB file')
        except:
                quitProgram(self, "No PDB File!")
                
        if self.pdbFile != None:
                cmd.load(self.pdbFile.name, 'PDB')
                print "Opening PDB file...", self.pdbFile.name
                get_PDBChains(self)
                for ch in self.chains:
                    print "Center of Geometry"
                    get_COG(self,ch)
                    draw_COG(self,ch)
                    calculate_tensor(self, ch)
                cmd.zoom()
        else:
                quitProgram(self,"No PDB file!")

def quitProgram(self,tit):
        tkMessageBox.showinfo(tit, "Quitting now...!")
        
def draw_COG(self, ch):
        com_object = 'COM_Chain_' + ch
        cmd.pseudoatom(object=com_object,pos=[self.cogx, self.cogy, self.cogz])
        cmd.show_as('spheres', com_object)
        cmd.color('yellow', com_object)

def get_PDBChains(self):
        self.chains = cmd.get_chains('PDB')
        for ch in self.chains:
                print "PDB has chain ", ch
        
def calculate_tensor(self, chain):
        selection = 'chain ' + chain
        object_name = 'GyrationAxes_' + chain
        model = cmd.get_model(selection, 1)
        nAtom = len(model.atom)
        
        xx,xy,xz,yy,yz,zz = 0,0,0,0,0,0
        x1,x2,x3,y1,y2,y3,z1,z2,z3 = 0,0,0,0,0,0,0,0,0
        x1n,x2n,x3n,y1n,y2n,y3n,z1n,z2n,z3n = 0,0,0,0,0,0,0,0,0
        x1a,x2a,x3a,y1a,y2a,y3a,z1a,z2a,z3a = 0,0,0,0,0,0,0,0,0
        x1an,x2an,x3an,y1an,y2an,y3an,z1an,z2an,z3an = 0,0,0,0,0,0,0,0,0
        
        self.evals, self.evecs, self.magnitude = [], [], []
        mcogx, mcogy, mcogz = 0,0,0
        
        for myAtom in model.atom:
                xx += (myAtom.coord[0] - self.cogx) * (myAtom.coord[0] - self.cogx)
                xy += (myAtom.coord[0] - self.cogx) * (myAtom.coord[1] - self.cogy)
                xz += (myAtom.coord[0] - self.cogx) * (myAtom.coord[2] - self.cogz)
                yy += (myAtom.coord[1] - self.cogy) * (myAtom.coord[1] - self.cogy)
                yz += (myAtom.coord[1] - self.cogy) * (myAtom.coord[2] - self.cogz)
                zz += (myAtom.coord[2] - self.cogz) * (myAtom.coord[2] - self.cogz)
                mcogx+= myAtom.coord[0] - self.cogx
                mcogy+= myAtom.coord[1] - self.cogy
                mcogz+= myAtom.coord[2] - self.cogz
        
        mcogx = mcogx/nAtom
        mcogy = mcogy/nAtom
        mcogz = mcogz/nAtom
        print "NEW COG : ", mcogx, mcogy, mcogz
        
        xx = xx/nAtom
        xy = xy/nAtom
        xz = xz/nAtom
        yy = yy/nAtom
        yz = yz/nAtom
        zz = zz/nAtom
        
        gyration_tensor = numpy.array([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]]);
        self.evals, self.evecs = numpy.linalg.eig(gyration_tensor)
        self.magnitude = numpy.sqrt(self.evals)
        
        print "\nMagnitude of Tensors :\n", self.magnitude, "\nEigen Vectors :\n", self.evecs
        
        x1,y1,z1 = self.magnitude[0] * self.evecs[0,0] + self.cogx, self.magnitude[0] * self.evecs[1,0] + self.cogy, self.magnitude[0] * self.evecs[2,0] + self.cogz
        x2,y2,z2 = self.magnitude[1] * self.evecs[0,1] + self.cogx, self.magnitude[1] * self.evecs[1,1] + self.cogy, self.magnitude[1] * self.evecs[2,1] + self.cogz
        x3,y3,z3 = self.magnitude[2] * self.evecs[0,2] + self.cogx, self.magnitude[2] * self.evecs[1,2] + self.cogy, self.magnitude[2] * self.evecs[2,2] + self.cogz
        
        x1n,y1n,z1n = self.magnitude[0] * -self.evecs[0,0] + self.cogx, self.magnitude[0] * -self.evecs[1,0] + self.cogy, self.magnitude[0] * -self.evecs[2,0] + self.cogz
        x2n,y2n,z2n = self.magnitude[1] * -self.evecs[0,1] + self.cogx, self.magnitude[1] * -self.evecs[1,1] + self.cogy, self.magnitude[1] * -self.evecs[2,1] + self.cogz
        x3n,y3n,z3n = self.magnitude[2] * -self.evecs[0,2] + self.cogx, self.magnitude[2] * -self.evecs[1,2] + self.cogy, self.magnitude[2] * -self.evecs[2,2] + self.cogz
        
        x1a,y1a,z1a = (self.magnitude[0] - 1) * self.evecs[0,0] + self.cogx, (self.magnitude[0] - 1) * self.evecs[1,0] + self.cogy, (self.magnitude[0] - 1) * self.evecs[2,0] + self.cogz
        x2a,y2a,z2a = (self.magnitude[1] - 1) * self.evecs[0,1] + self.cogx, (self.magnitude[1] - 1) * self.evecs[1,1] + self.cogy, (self.magnitude[1] - 1) * self.evecs[2,1] + self.cogz
        x3a,y3a,z3a = (self.magnitude[2] - 1) * self.evecs[0,2] + self.cogx, (self.magnitude[2] - 1) * self.evecs[1,2] + self.cogy, (self.magnitude[2] - 1) * self.evecs[2,2] + self.cogz
        
        x1an,y1an,z1an = (self.magnitude[0] - 1) * -self.evecs[0,0] + self.cogx, (self.magnitude[0] - 1) * -self.evecs[1,0] + self.cogy, (self.magnitude[0] - 1) * -self.evecs[2,0] + self.cogz
        x2an,y2an,z2an = (self.magnitude[1] - 1) * -self.evecs[0,1] + self.cogx, (self.magnitude[1] - 1) * -self.evecs[1,1] + self.cogy, (self.magnitude[1] - 1) * -self.evecs[2,1] + self.cogz
        x3an,y3an,z3an = (self.magnitude[2] - 1) * -self.evecs[0,2] + self.cogx, (self.magnitude[2] - 1) * -self.evecs[1,2] + self.cogy, (self.magnitude[2] - 1) * -self.evecs[2,2] + self.cogz
        
        self.obj = []
        self.obj.extend([ 9.0, self.cogx, self.cogy, self.cogz, x1a, y1a, z1a, 0.25, 1, 0, 0, 1, 0, 0 ])
        self.obj.extend([ 9.0, self.cogx, self.cogy, self.cogz, x2a, y2a, z2a, 0.25, 0, 1, 0, 0, 1, 0 ])
        self.obj.extend([ 9.0, self.cogx, self.cogy, self.cogz, x3a, y3a, z3a, 0.25, 0, 0, 1, 0, 0, 1 ])
        self.obj.extend([CONE, x1a, y1a, z1a, x1, y1, z1, 0.5, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0])
        self.obj.extend([CONE, x2a, y2a, z2a, x2, y2, z2, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0])
        self.obj.extend([CONE, x3a, y3a, z3a, x3, y3, z3, 0.5, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        self.obj.extend([ 9.0, self.cogx, self.cogy, self.cogz, x1an, y1an, z1an, 0.25, 1, 0, 0, 1, 0, 0 ])
        self.obj.extend([ 9.0, self.cogx, self.cogy, self.cogz, x2an, y2an, z2an, 0.25, 0, 1, 0, 0, 1, 0 ])
        self.obj.extend([ 9.0, self.cogx, self.cogy, self.cogz, x3an, y3an, z3an, 0.25, 0, 0, 1, 0, 0, 1 ])
        self.obj.extend([CONE, x1an, y1an, z1an, x1n, y1n, z1n, 0.5, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0])
        self.obj.extend([CONE, x2an, y2an, z2an, x2n, y2n, z2n, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0])
        self.obj.extend([CONE, x3an, y3an, z3an, x3n, y3n, z3n, 0.5, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        cmd.load_cgo(self.obj,object_name)
        
def get_COG(self, chain):
        selection = 'chain ' + chain
        model = cmd.get_model(selection,1)
        self.cogx,self.cogy,self.cogz=0,0,0
        nAtom = len(model.atom)
        
        for myAtom in model.atom:
                self.cogx+= myAtom.coord[0]
                self.cogy+= myAtom.coord[1]
                self.cogz+= myAtom.coord[2]
       
        self.cogx=self.cogx/nAtom
        self.cogy=self.cogy/nAtom
        self.cogz=self.cogz/nAtom
        print "X: ", self.cogx, "Y: ", self.cogy, "Z: ", self.cogz
        