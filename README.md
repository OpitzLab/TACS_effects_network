# TACS_effects_network
We investigated the effects of tACS intensities and frequencies on computational neuronal network model of two-compartment neurons.

## Information about the cell model used

## Software requirements and instructions
The data were generated using NEURON Environment 8.0.0. You can download this software here https://neuron.yale.edu/neuron/. 
Assuming NEURON has been installed, the MOD files in the 'mod/' folder need to be compiled using mknrndll. Copy then nrnmech.dll from 'mod/' up one level to the main folder.

NetPyNE is required to run the network model. It can be downloaded by following http://www.netpyne.org/install.html.
If you successfully installed MPI (e.g. OpenMPI) and NEURON with MPI support, you can simulate the model in parallel using multiple cores/processors by typing:
mpiexec -n 4 -python run_2C_tACS_network_par.py

## License
The code is licensed with the CC-BY-NC-SA license.

## Correspondence
aopitz -at- umn.edu (Alex Opitz)
