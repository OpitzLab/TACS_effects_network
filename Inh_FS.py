# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 13:45:27 2021

TEMPLATE FILE FOR FAST-SPIKING CORTICAL INTERNEURON
one compartment model and currents derived from:

   Pospischil, M., Toledo-Rodriguez, M., Monier, C., Piwkowska, Z., 
   Bal, T., Fregnac, Y., Markram, H. and Destexhe, A.
   Minimal Hodgkin-Huxley type models for different classes of
   cortical and thalamic neurons.
   Biological Cybernetics 99: 427-441, 2008.

Two compartment model modified by Zhihe Zhao
	- two compartment model
	- passive
	- HH: Traub
    
@author: zhaoz
"""

from neuron import h
from math import pi

class Inh:
  def __init__ (self,alpha=2.0):
    self.soma = soma = h.Section(name='soma',cell=self)
    self.dend = dend = h.Section(name='dend',cell=self)
    self.dend.connect(self.soma,0.5,0) #   connect dend(0), soma(0.5)
    self.alpha = alpha
    for sec in [self.soma,self.dend]:
        sec.insert('pas')
        sec.e_pas = -70
        sec.g_pas = 0.00015*(.5) 
        sec.insert('extracellular') # extracellular mech
        sec.insert('xtra') # xtra.mod
    self.initsoma()
    self.initdend()
    
  def initsoma (self):
    soma = self.soma
    soma.nseg = 1
    soma.diam = 50/self.alpha
    soma.L = 50/self.alpha
    soma.cm = 0.0005*(self.alpha)**2
    soma.Ra = 35.4
    soma.insert('hh2') # Hodgin-Huxley INa and IK 
    soma.vtraub_hh2 = -55 # resting Vm
    soma.ek = -90 # K reversal potential(mV)
    soma.ena = 50 # Na reversal potential(mV)
    soma.gnabar_hh2 = 0.05*(self.alpha**2)
    soma.gkbar_hh2 = 0.01*(self.alpha**2) #spike duration of interneuron 

  def initdend (self):
    dend = self.dend
    dend.nseg = 1
    dend.diam = 10.0/pi 
    dend.cm = 1
    dend.L = 200 # microns
    dend.Ra = 100