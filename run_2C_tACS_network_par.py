# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 11:02:31 2022

@author: zhaoz
"""

# from IPython import get_ipython
# get_ipython().magic('reset -f') 
from mpi4py import MPI
import numpy as np
from netpyne import specs, sim
from neuron import h # Import NEURON
from math import pi

from NetWeight import *

pc = h.ParallelContext()
pc.timeout(0)

nettype = 1 #1 alpha oscillation; 2 low beta oscillation (~15 Hz); 0 no clear oscillation 
netweight = NetworkWeight(nettype)
# STIM_TYPE = 1
twait = 500 # ms ; wait time for the network to reach stable state
simdur = 240000 # ms ; simulation duration(ms)
DEL = twait+0.5*simdur    #stim delay (ms)
DUR = 0.5*simdur          #stim duration (ms)
AMP = 0                   #stim amplitude (mV/mm)
offset = 0
FREQ = 10                 #stim frequency (Hz)
dt = 0.05 #ms

# defines the direction of the E-field (y )
Ex = 0
Ey = -1
Ez = 0
E1 = [Ex,Ey,Ez] # E-field strength mV/mm

a_PY = 5 # the scaling factor of E-field coupling of PY
a_IN = 1.5 # the scaling factor of E-field coupling of IN


# h.load_file("stimWaveform.hoc")

# Network parameters
netParams = specs.NetParams()  # object of class NetParams to store the network parameters
netParams.sizeX = 200 # x-dimension (horizontal length) size in um
netParams.sizeY = 1000 # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = 20 # z-dimension (horizontal length) size in um

exc_rho = 200 
exc_soma_L = 30 #um
exc_soma_diam = 20
exc_dend_L = exc_rho * exc_soma_L/3
inh_soma_L = 25
inh_soma_diam = 25
inh_dend_L = 200
inh_dend_diam = 10/pi

cellRule = netParams.importCellParams(
        label='PYR_Mainen_rule',
        conds={'cellType': 'PYR', 'cellModel': 'Mainen'},
        fileName='mainen.py',
        cellName='PYR2')
cellRule['secs']['soma']['vinit'] = -70 #initial membrane voltage (mV)
cellRule['secs']['soma']['geom']['pt3d'] = []
cellRule['secs']['soma']['geom']['pt3d'].append((0, 0, 0, 20.0))
cellRule['secs']['soma']['geom']['pt3d'].append((0, 30, 0, 20.0))

cellRule['secs']['dend']['geom']['pt3d'] = []
cellRule['secs']['dend']['geom']['pt3d'].append((0, 30, 0, 10.0/pi))
cellRule['secs']['dend']['geom']['pt3d'].append((0, 200*10+30, 0, 10.0/pi))
# netParams.cellParams['PYR_rule'] = cellRule

cellRule = netParams.importCellParams(
        label='Inh_rule',
        conds={'cellType': 'Inh', 'cellModel': 'Inh_FS'},
        fileName='Inh_FS.py',
        cellName='Inh')
cellRule['secs']['soma']['geom']['pt3d'] = []
cellRule['secs']['soma']['geom']['pt3d'].append((0, 40, 0, 25))
cellRule['secs']['soma']['geom']['pt3d'].append((0, 65, 0, 25))
cellRule['secs']['dend']['geom']['pt3d'] = []
cellRule['secs']['dend']['geom']['pt3d'].append((0, 65, 0, 10.0/pi))
cellRule['secs']['dend']['geom']['pt3d'].append((0, 265, 0, 10.0/pi))
# netParams.cellParams['INT_rule'] = cellRule

# Population parameters
num_exc = 800
num_inh = 200
# 
netParams.popParams['exc'] = {'cellType': 'PYR', 'numCells': num_exc,'xRange': [100, 1100],'yRange': [600, 1000],'zRange':[-500, 500],'cellModel': 'Mainen'}
netParams.popParams['Inh'] = {'cellType': 'Inh', 'numCells': num_inh, 'xRange': [100, 1100],'yRange': [700, 800],'zRange':[-500, 500],'cellModel': 'Inh_FS'}

# Synaptic mechanism parameters
netParams.synMechParams['exc'] = {'mod': 'Exp2Syn', 'tau1': 2, 'tau2': 10, 'e': 0}  # excitatory synaptic mechanism

#Synaptic mechanism parameters are based on C.Rusu 2014 Brain Stim
netParams.synMechParams['AMPA'] = {'mod': 'Exp2Syn', 'tau1': 0.2, 'tau2': 1.7, 'e': 0}  # excitatory synaptic mechanism
netParams.synMechParams['NMDA'] = {'mod': 'Exp2Syn', 'tau1': 2, 'tau2': 26, 'e': 0}  # excitatory synaptic mechanism
netParams.synMechParams['GABAA'] = {'mod': 'Exp2Syn', 'tau1': 0.3, 'tau2': 2.5, 'e': -70}  # excitatory synaptic mechanism
netParams.synMechParams['GABAB'] = {'mod': 'Exp2Syn', 'tau1': 45.2, 'tau2': 175.16, 'e': -90}  # excitatory synaptic mechanism

# Cell connectivity rules


netParams.connParams['exc->exc'] = { #  S -> M label
        'preConds': {'pop': 'exc'},   # conditions of presyn cells
        'postConds': {'pop': 'exc'},  # conditions of postsyn cells
        'sec' : 'dend',         # connection to dendrite
        'probability': 0.1,         # probability of connection
        'weight': [netweight['E2E_AMPA'], netweight['E2E_NMDA']],             # synaptic weight
        'delay': [1,1],                 # transmission delay (ms)
        'synMech': ['AMPA','NMDA']}           # synaptic mechanism

netParams.connParams['exc->Inh'] = { #  S -> M label
        'preConds': {'pop': 'exc'},   # conditions of presyn cells
        'postConds': {'pop': 'Inh'},  # conditions of postsyn cells
        'sec' : 'dend',         # connection to dendrite
        'probability': 0.1,         # probability of connection
        'weight': [netweight['E2I_AMPA'], netweight['E2I_NMDA']],             # synaptic weight
        'delay': [1,1],                 # transmission delay (ms)
        'synMech': ['AMPA','NMDA']}           # synaptic mechanism

netParams.connParams['Inh->exc'] = { #  S -> M label
        'preConds': {'pop': 'Inh'},   # conditions of presyn cells
        'postConds': {'pop': 'exc'},  # conditions of postsyn cells
        'sec' : 'soma', # Inh connection to soma
        'probability': 0.15,         # probability of connection
        'weight': netweight['I2E_GABA'],          # synaptic weight
        'delay': 1,                 # transmission delay (ms)
        'synMech': 'GABAA'}           # synaptic mechanism

netParams.connParams['Inh->Inh'] = { #  S -> M label
        'preConds': {'pop': 'Inh'},   # conditions of presyn cells
        'postConds': {'pop': 'Inh'},  # conditions of postsyn cells
        'sec' : 'soma', # Inh connection to soma
        'probability': 0.1,         # probability of connection
        'weight':netweight['I2I_GABA'],             # synaptic weight
        'delay': 1,                 # transmission delay (ms)
        # 'synMech': ['GABAA','GABAB']
        'synMech': 'GABAA'
        }           # synaptic mechanism

# Stimulation parameters
netParams.stimSourceParams['bkg'] = {'type': 'NetStim', 'interval' : 1000/50, 'noise': 1}
netParams.stimTargetParams['bkg->PYR'] = {'source': 'bkg', 'conds': {'cellType': 'PYR'}, 'sec':'dend', 'weight': netweight['PoissonE'], 'delay': 0, 'synMech': 'exc'}
netParams.stimTargetParams['bkg->Inh'] = {'source': 'bkg', 'conds': {'cellType': 'Inh'}, 'sec':'dend', 'weight': netweight['PoissonI'], 'delay': 0, 'synMech': 'exc'}

# parameters for IN orientation
np.random.seed(0) #set seed for rand.uniform
theta = np.random.uniform(0,pi,size=num_inh)
phi = np.random.uniform(-pi,pi,size=num_inh)

# Simulation options
simConfig = specs.SimConfig()       # object of class SimConfig to store simulation configuration

simConfig.duration = twait+simdur         # Duration of the simulation, in ms
simConfig.dt = dt               # Internal integration timestep to use
simConfig.hParams = {'celsius': 37}
simConfig.seeds['stim'] = 1
simConfig.verbose = False           # Show detailed messages 
simConfig.saveCellConns = True
simConfig.recordTraces = {'V_soma':{'sec':'soma','loc':0.5,'var':'v'}}  # Dict with traces to record
simConfig.recordTraces['tACS'] = {'sec':'dend','loc':0.5,'var':'ex_xtra'} # Dict with traces to record
# simConfig.recordTraces['i_AMPA'] = {'sec': 'dend', 'loc': 0.5,'synMech': 'AMPA','var': "i"} # record the potassium current at the proximal end of the first 'branch' section (branch[0])
# simConfig.recordTraces['i_GABAA'] = {'sec': 'soma', 'loc': 0.5,'synMech': 'GABAA','var': "i"} # record the potassium current at the proximal end of the first 'branch' section (branch[0])
simConfig.recordStep = 0.05            # Step size in ms to save data (eg. V traces, LFP, etc)
simConfig.saveDataInclude = ['simData','netParams'] 
simConfig.filename = 'tACS_network_results'         # Set file output name
simConfig.savePickle = True       # Save params, network and sim output to pickle file
simConfig.recordLFP = [[x, y, 0] for y in range(200, 1000, 200) for x in [600]]


simConfig.analysis['plotRaster'] = {'saveFig': True,'timeRange': [twait,simConfig.duration]}                  # Plot a raster
simConfig.analysis['plotTraces'] = {'include': [0, 91], 'saveFig': True,'overlay': False, 'oneFigPer': 'trace','showFig': True}  # Plot recorded traces for this list of cells
# simConfig.analysis['plot2Dnet'] = {'saveFig': False}                   # plot 2D cell positions and connections
simConfig.analysis['plotLFP'] = {'includeAxon': False, 'figSize': (6,10), 'timeRange': [twait,simConfig.duration], 'saveFig': True,'filtFreq':300, 'plots': ['timeSeries','PSD','location']} 
# simConfig.analysis['plotSpikeStats'] = {'timeRange': [twait,simConfig.duration], 'graphType': 'histogram', 'bins': 25, 'stats': 'rate' , 'saveFig' : True }

# Create network and run simulation
sim.initialize(                     # create network object and set cfg and net params
        simConfig = simConfig,          # pass simulation config and network params as arguments
        netParams = netParams)
sim.net.createPops()                # instantiate network populations
sim.net.createCells()               # instantiate network cells based on defined populations
sim.net.connectCells()              # create connections between cells based on params
sim.net.addStims()                  # add stimulation
cellsOnCurrentRank = sim.net.cells
# print("I'm rank {sim.rank} and I have {len(cellsOnCurrentRank)} cells")

#tACS parameters
h.load_file('stdrun.hoc')

h('{DEL = 300*1000 DUR = 600*1000 AMP= 1 offset = 0 FREQ = 10 STIM_TYPE = 1}')
h.DEL = DEL
h.DUR = DUR
h.AMP = AMP
h.offset = offset
h.FREQ = FREQ
h.STIM_TYPE = 1
h.tstop = twait + simdur
h.dt = dt

h.load_file("stimWaveform.hoc")

#randomize interneuron orientation in the Network
for n in range (len(cellsOnCurrentRank)):
    if 'Inh' in cellsOnCurrentRank[n].tags['cellType']:
        sec_dend = cellsOnCurrentRank[n].secs['dend']['hObj']
        phi1 = phi[cellsOnCurrentRank[n].gid-num_exc]
        theta1 = theta[cellsOnCurrentRank[n].gid-num_exc]
        sec_dend.pt3dchange(1,sec_dend.x3d(1)+inh_dend_L*np.cos(phi1)*np.sin(theta1),
                            sec_dend.y3d(1)-inh_dend_L+inh_dend_L*np.cos(theta1),
                            sec_dend.z3d(1)+inh_dend_L*np.sin(theta1)*np.sin(phi1),
                            sec_dend.diam3d(1))

#setpointer to link e_extracellular to ex_xtra
for sec in h.allsec():
    h.ion_style("ca_ion",0,1,0,0,0, sec = sec) # using an ohmic current rather than GHK equation
    if h.ismembrane('xtra', sec=sec) and h.ismembrane('extracellular', sec=sec):
        # print('sec has xtra and extracellular')
        h.setpointer(sec(0.5)._ref_e_extracellular, 'ex', sec(0.5).xtra)
        
#Assigne 3d coordinates to xtra       
for m in range (len(cellsOnCurrentRank)):
    secdend = sim.net.cells[m].secs['dend']['hObj']
    secsoma = sim.net.cells[m].secs['soma']['hObj']

    secdend.x_xtra = secdend.x3d(1)
    secdend.y_xtra = secdend.y3d(1)
    secdend.z_xtra = secdend.z3d(1)
    
    secsoma.x_xtra = secsoma.x3d(0)
    secsoma.y_xtra = secsoma.y3d(0)
    secsoma.z_xtra = secsoma.z3d(0)


# #Apply tACS to the network
for j in range (len(cellsOnCurrentRank)):
    if 'PYR' in cellsOnCurrentRank[j].tags['cellType']:
        cell=cellsOnCurrentRank[j]
        sec = cell.secs['dend']['hObj']
        if h.ismembrane('xtra',sec=sec):
            sec.es_xtra = -(Ex*sec.x_xtra + Ey*sec.y_xtra + Ez*sec.z_xtra)*1e-3/a_PY
        sec = cell.secs['soma']['hObj']
        if h.ismembrane('xtra',sec=sec):
            sec.es_xtra = -(Ex*sec.x_xtra + Ey*sec.y_xtra + Ez*sec.z_xtra)*1e-3/a_PY
    elif 'Inh' in cellsOnCurrentRank[j].tags['cellType']:
        cell=cellsOnCurrentRank[j]
        sec = cell.secs['dend']['hObj']
        if h.ismembrane('xtra',sec=sec):
            sec.es_xtra = -(Ex*sec.x_xtra + Ey*sec.y_xtra + Ez*sec.z_xtra)*1e-3/a_IN
        sec = cell.secs['soma']['hObj']
        if h.ismembrane('xtra',sec=sec):
            sec.es_xtra = -(Ex*sec.x_xtra + Ey*sec.y_xtra + Ez*sec.z_xtra)*1e-3/a_IN     
        
        
# set tACS stimulation         
h('setstim(DEL,DUR,AMP,FREQ,offset)')

sim.setupRecording()                # setup variables to record for each cell (spikes, V traces, etc)
sim.runSim()                        # run parallel Neuron simulation
sim.gatherData()                    # gather spiking data and cell info from each node
sim.saveData()                      # save params, cell info and sim output to file (pickle,mat,txt,etc)
# sim.runSimWithIntervalFunc(100, sim.intervalSave)
# sim.gatherData()
# sim.intervalSimulate(100)
sim.analysis.plotData()             # plot spike raster   

# #plot spike stats
sim.analysis.plotSpikeStats(include=['exc','Inh'], statDataIn={}, timeRange=[twait,simConfig.duration], graphType='histogram', stats=['rate'],bins = 25,legendLabels= ['Inh','exc'], saveFig = True)

# import matplotlib.pyplot as plt
# plt.close('all')
# import scipy.signal
# from scipy.stats import kurtosis

# fs = 1000/simConfig.dt
# netLFP = sim.simData['LFP'][np.int(twait*fs/1000):np.int(simConfig.duration*fs/1000),3]
# # netLFP = sim.simData['LFP'][:,7]

# (freq, netPSD)= scipy.signal.welch(netLFP,  fs, nperseg= fs)
# maxnetPSDFreq = float(freq[np.where(netPSD == np.amax(netPSD))])

# plt.figure
# plt.plot(freq, netPSD)
# plt.xlim([0, 100])
# plt.xlabel('frequency [Hz]')
# plt.ylabel('PSD [V**2/Hz]')
# plt.savefig("tACS_network_PSD.png")
# plt.close()

