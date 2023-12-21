# Binding_Polynomial_Analysis
## Basic overview and theory
Codes from recent paper to conduct analysis on standard MD simulations, and fit ligand binding data to a polynomial model.
This contains two directories:

**01_generate_histograms** contains a sample trajectory and topology file from the paper, an MD simulation of a Fab fragment of MEDI-4212 and 25 free lysine amino acids. 
The anlysis script, get_distmatrix_nbound_vs_time.py calculates all contacts between each protein residue and each free lysine molecule at every step of the trajectory. 
To use the script on your own simulation trajectories, you will need to modify the script with the appropriate file names, as well as the appropriate amino acid residue indices which you want to calculate binding to. 
The main output from this directory will be a histogram of the number of lysine molecules bound to the given residues on the protein surface over the time of the simulation.

**02_fit_histograms** contains a script to fit the binomial theory developed in this work which we detail below.

We derive a probability distribution of number of bound excipient molecules from first principles, and use this theory to fit the data. The fraction of protein molecules, with _m_ bound excipient molecules is expressed as $p(PE_m)$. From first principles, we find this to be a binomial distribution, having the functional form:

### $p(PE_m)=\frac{{N \choose m} (K_e c)^m}{(1+K_e c)^N}$

where $N$ is the maximum number of excipients that can bind to the interface, $K_e$ is the association constant of a single excipient molecule with a site on the protein surface, and c is the concentration of excipient in solution. 
Within this study, we find that N and $K_e$ usually do not have unique solutions, and are in fact inversely related to each other.
Thus, we define a proxy parameter: $a=N \times K_e$, which does have a unique solution, and is in fact a good proxy variable to describe how effective an excipient molecule is at destabilizing protein-protein interactions.
The equation that we fit to in the script is as follows:

### $p(PE_m)=\frac{{N \choose m} (\frac{a}{N} c)^m}{(1+\frac{a}{N} c)^N}$

and optimize the fit to the $a$ and $N$ parameters.

## To Execute
1. Ensure your python has [MDTraj](https://www.mdtraj.org/1.9.8.dev0/index.html), [Numpy](https://numpy.org/) and [Scipy](https://scipy.org/) installed. The script also uses [MatPlotLib](https://matplotlib.org/) to generate some plots.

2. Download or clone the repository

3. Navigate to the first directory
```console
cd Binding_Polynomial_Analysis/01_generate_histograms
```

4. Execute the first script
```console
python get_distmatrix_nbound_vs_time.py
```

5. If the script runs correctly, it should display two plots, one heatmap, and one two-panel graph showing the timeseries of number of bound excipients, and the histogram.
There should also be three output files, one .npy and two .dat files. You can also uncomment two lines in the code to make it output .csv files to plot in Excel or Numbers.
To run this with your own trajectories, change the file names in the script, as well as the list of residues identified to be in the PPI interface, and the chain IDs to agree with your system.

6. Once you have generated histograms, you can then copy them to the second directory
```console
cp *histogram.dat ../02_fit_histograms
```
You will want to rename them or modify the script in that directory so that they are recognized by the script.

7. Once the histogram files are named appropriately, and the script is set up to recognize them, execute the script in this directory:
```console
python fit_binomial_N_a.py
```
This will conduct the optimization, and plot the data for each concentration you have tested, as well as displaying the optimized $a$ and $N$ parameters.
