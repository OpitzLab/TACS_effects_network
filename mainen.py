# translated from /u/samn/npredict/geom_mainen.hoc // $Id: geom_mainen.hoc,v 1.3 2007/04/16 14:42:33 samn Exp $
from neuron import h
from math import pi

class PYR2:
  def __init__ (self,ID=0,ty=0,col=0,rho=200.0,kappa=10.0,soma_pas=True):
    self.ID=ID
    self.ty=ty
    self.col=col
    self.soma_pas = soma_pas
    self.soma = soma = h.Section(name='soma',cell=self)
    self.dend = dend = h.Section(name='dend',cell=self)
    self.dend.connect(self.soma,0.5,0) #   connect dend(0), soma(0.5)
    self.rho = rho # dendritic to axo-somatic area ratio
    self.kappa = kappa # coupling resistance (Mohm)
    for sec in [self.soma,self.dend]:
        sec.insert('k_ion')
        sec.insert('na_ion')
        sec.insert('ca_ion')
        sec.ek = -90 # K+ current reversal potential (mV)
        sec.ena = 60 # Na+ current reversal potential (mV)
        sec.eca = 140 # Ca2+ current reversal potential (mV)
        h.ion_style("ca_ion",0,1,0,0,0) # using an ohmic current rather than GHK equation
        sec.Ra= 100
        sec.insert('extracellular') # extracellular mech
        sec.insert('xtra') # xtra.mod
    self.initsoma()
    self.initdend()

  def initsoma (self,alpha1 = 2*pi,alpha = 3):
    self.alpha = alpha #scaling factor for morphology and properties 
    soma = self.soma
    soma.nseg = 1
    soma.diam = 10.0/pi*alpha1
    soma.L = 10*alpha
    soma.cm = 0.75*6/(alpha*alpha1)
    soma.insert('na') # na.mod
    soma.insert('kv') # kv.mod
    soma.gbar_na = (30e3)/(alpha*alpha1*0.82)
    soma.gbar_kv = (1.5e3)/(alpha*alpha1*1)
    # soma.insert('extracellular') # extracellular mech
    # soma.insert('xtra') # xtra.mod
    if self.soma_pas:
      soma.insert('pas')
      soma.e_pas=-70
      soma.g_pas=1/3e4 
    # soma.pt3dadd(0,0,0,soma.diam)
    # soma.pt3dadd(0,10,0,soma.diam)

  def initdend (self):
    dend = self.dend
    dend.nseg = 1
    dend.diam = 10.0/(2*pi)
    self.config()
    dend.cm = 0.75
    dend.insert('na') # na.mod
    dend.insert('km') # km.mod
    dend.insert('kca') # kca.mod
    dend.insert('ca') # ca.mod
    dend.insert('cad') # cad.mod
    dend.insert('pas')
    # dend.insert('extracellular') # extracellular mech
    # dend.insert('xtra') # xtra.mod
    dend.eca = 140
    h.ion_style("ca_ion",0,1,0,0,0) # already called before
    dend.e_pas = -70 # only dendrite has leak conductance - why?
    dend.g_pas = 1/3e4 # only dendrite has leak conductance
    dend.gbar_na= 15
    dend.gbar_ca = 0.3 # high voltage-activated Ca^2+
    dend.gbar_km = 0.1 # slow voltage-dependent non-inactivating K+
    dend.gbar_kca = 3  # slow Ca^2+-activated K+

        
        
  def config (self):
    self.dend.L  = self.rho*self.soma.L/self.alpha # dend area is axon area multiplied by rho
    self.dend.Ra = self.dend.Ra*self.kappa/self.dend(0.5).ri() # axial resistivity is adjusted to achieve

  # resets cell to default values
  def todefault(self):
    self.rho = 100
    self.kappa = 10
    self.config()
    self.soma.gbar_na = 30e3
    self.soma.gbar_kv = 1.5e3
    self.dend.g_pas = 1/3e4
    self.dend.gbar_na = 15
    self.dend.gbar_ca = 0.3
    self.dend.gbar_km = 0.1
    self.dend.gbar_kca = 3
    self.soma.ek = dend.ek = -90
    self.soma.ena = dend.ena = 60
    self.soma.eca = dend.eca = 140

