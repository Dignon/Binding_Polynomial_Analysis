#!/usr/bin/python
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

"""This script takes an MD trajectory and a topology file and outputs a distance matrix
showing the distance between every excipient molecule with every residue of the 
protein at every timestep. Thus the output is a 3D numpy array with dimensions
nframes x nresidues x nexcipients. The data is visualized at the end by averaging
over the time axis (0), showing the average distance between each excipient molecule
and each amino acid residue through the whole simulation. Distance is defined as
the minimum distance between heavy atoms in the excipient molecule and the amino
acid residue being looked at. Depending on the length and size of your simulation,
this might take a while to run. You can run a quick preliminary test by increasing
the stride of the analysis, e.g. to 100, to ensure the script is working properly.
After calculating the distance matrix, it then takes the identified residues from
the PPI interface, and determines how may excipient molecules are bound to this 
region at each timestep in the simulation. Finally, it generates a histogram of
these numbers which can then be used to fit to the binomial equation from the
source article.

The output files generaed by this script are:
distance_matrix.npy: A 3D numpy array, dimension nframes x nresidues x nexcipients, containing the minimum heavy atom distance between each amino acid residue in the protein, and each excipient molecule in the system at every time step.
n_interface_contacts.{dat,csv}: Timeseries of number of excipient molecules bound to the PPI interface at every frame of the simulation.
contact_histogram.{dat,csv}: A histogram of the timeseries in previously mentioned file.

Please cite:
Dignon GL, Dill KA, bioRxiv 2023 DOI:
"""

## Change these parameters to include your protein file, as well as indices
# you would like to consider PPI regions
Trajectory_file = "Traj_nowater_short.dcd" #This could be xtc, dcd, etc. or a list of files
Topology_file = "updated_nowater.pdb" #Make sure this has chain identifiers for protein and excipient molecules
start_frame = 0     # First frame of the trajectory to read
end_frame = 5000    # Last frame of the trajectory to read
stride = 1          # Only analyze every <stride> frames in the trajectory 
cutoff_dist = 0.8 #nm,  cutoff distance to consider a protein-excipient contact

## Determine which residues in your protein are involved in PPI and number them
# starting at 0. Change the indices list to your list of residues
PPI_indices = [0, 2, 4, 22, 25, 50, 53, 54, 67, 68, 69, 70, 270, 273, 319] # residues in example Fab protein, from model 3 with 4A^2 SASA cutoff definition for PPI residues

## Define atom indices in different groups
traj = md.load(Trajectory_file, top=Topology_file)
protein_idx = traj.topology.select("chainid 0 to 1") # chains A and B are the Fab. If your protein has a different number of chain segments (not 2) make sure to change this line, and the following line as well.
excipient_idx = traj.topology.select("chainid 2") # chain C is the excipient molecules

## Make partial trajectory objects
protein = traj.atom_slice(protein_idx, inplace=False)
excipient = traj.atom_slice(excipient_idx, inplace=False)

## Get list of pairs that will be compared to get protein-excipient interactions
n_prot_res = protein.topology.n_residues
n_excipient_res = excipient.topology.n_residues

# protein column [1,1,1...1,2,2...,n_res,n_res]
prot_col = np.array([[i]*n_excipient_res for i in range(n_prot_res)])
prot_col = prot_col.flatten()

# excipient column [1,2,3,...n_exc,1,2,3,...,n_exc]
exc_col = list(range(n_prot_res, n_prot_res + n_excipient_res)) * n_prot_res

# set up pair list and copy in columns
pair_list = np.zeros((n_prot_res*n_excipient_res, 2))
pair_list[:,0] = prot_col
pair_list[:,1] = exc_col
print(pair_list)

## Calculate contacts in between excipient molecules and protein residues
contacts, indices = md.compute_contacts(traj[start_frame:end_frame:stride], contacts=pair_list, scheme='closest-heavy', periodic=True)

distance_map = contacts.reshape(len(traj[start_frame:end_frame:stride]), n_prot_res, n_excipient_res)
print(distance_map.shape)
np.save('distance_map.npy', distance_map)

im = plt.imshow(distance_map.mean(0), cmap='jet', aspect='auto')
cbar = plt.colorbar()
plt.xlabel('Excipient ID')
plt.ylabel('Residue ID')
cbar.ax.set_ylabel('Average Distance (nm)')
plt.show()

## Convert the distance map to a contacts map
contacts = distance_map <= cutoff_dist

## Isolate just the contacts occurring with the identified residues in the PPI interface
Interface_contacts = contacts[:,PPI_indices,:] # 3D array size: nframes x len(PPI_indices) x nexcipient

## Find how many unique excipient molecules are bound to the interface in each frame
n_contacts = (Interface_contacts.sum(1) > 0).sum(1) # 1D array size: nframes

## Plot the results
plt.subplot(1,2,1)
plt.plot(np.arange(len(n_contacts))*stride, n_contacts)
plt.xlabel('Frame')
plt.ylabel('Number of bound excipient molecules')
plt.ylim(-0.5,n_contacts.max()+0.5)

## Format results and output into text file
outp = np.zeros((distance_map.shape[0],2))
outp[:,0] = np.arange(distance_map.shape[0])
outp[:,1] = n_contacts
np.savetxt('n_interface_contacts.dat', outp, fmt='%-8i %-4i')
#np.savetxt('n_interface_contacts.csv', outp, fmt='%i', delimiter=',') # uncomment if you want csv file to import into Excel

## Histogram the contact data and output to file
bins, hist = np.unique(n_contacts, return_counts=True)
out_hist = np.zeros((len(bins),2))
out_hist[:,0] = bins
out_hist[:,1] = hist
np.savetxt('contact_histogram.dat', out_hist, fmt='%-4i %-8i')
#np.savetxt('contact_histogram.csv', out_hist, fmt='%i', delimiter=',')

## Plot the histogram alongside the timeseries
plt.subplot(1,2,2)
plt.barh(bins, hist)
plt.ylabel('# Bound excipients')
plt.xlabel('# Frames')
plt.ylim(-0.5,n_contacts.max()+0.5)
plt.show()

