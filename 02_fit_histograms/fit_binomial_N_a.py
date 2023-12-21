#!/usr/bin/python
import numpy as np
from scipy.special import comb
from scipy.optimize import curve_fit
import scipy.optimize as optimize
import matplotlib.pyplot as plt

"""This script will take data from simulations and fit it to a binomial. For the 
example here, we will show a multi-curve fit where we possess three histograms of
a single protein-exciipent combination at different excipient concentrations. These
will be fit simultaneously to the binomial theory developed in this work.
"""

## First define the functions that will be used to fit the data
def BindingBinomialNaRelation(n, cE, a, N):
    return comb(N, n) * ((a*cE/(N))**(n))

def NormedBindingBinomialNaRelation(param, refdat, cElist):
    ## Unpack parameters
    a = param[0]
    N = param[1]

    ## Ref should contain the list of data sets being fitted to (n, pPEn)
    loss = 0.
    for i, cE in enumerate(cElist):
        nlist = refdat[i][:,0]
        pPEnlist = refdat[i][:,1]
        unnormed = BindingBinomialNaRelation(nlist, cE, a, N)
        normed = unnormed / unnormed.sum()
        loss += np.sqrt(((normed - pPEnlist)**2).sum())
    return loss

## Set the list of excipient concentrations used to obtain the histograms
cElist = np.array([0.077,0.155,0.233]) 

## Import each data set
fig, ax = plt.subplots(dpi=150)
ref_data = []
for cE in cElist:
    ref_data.append(np.loadtxt('Lys_%03d_hist.dat'%(cE*1000)))

    ## Normalize the histograms
    ref_data[-1][:,1:] = ref_data[-1][:,1:] / ref_data[-1][:,1:].sum(0)
    ref_data[-1] = ref_data[-1][:,[0,1]]

    ## Get rid of non-zero entries
    ref_data[-1] = ref_data[-1][ref_data[-1][:,1]!=0]
print(ref_data)

## Set initial guesses for the a and N parameters, kE is mostly independent of these
param_ini = [40.0,25] # a, N

## Use optimizing protocol to fit multiple curves and minimize "loss" or deviation between data and binomial fit
op = optimize.basinhopping(NormedBindingBinomialNaRelation, param_ini, minimizer_kwargs={'args': (ref_data,cElist), 'method': "L-BFGS-B", 'bounds': ([0.1,100],[1,100])})
print(op)
a, N = op.x

## Plot the data and the resulting binomial fits and fitted parameters
clist = plt.cm.jet(np.linspace(0,1,len(cElist)))
for i, cE in enumerate(cElist):
    ax.scatter(ref_data[i][:,0], ref_data[i][:,1], c=clist[i])
    conc_PEn = BindingBinomialNaRelation(ref_data[i][:,0], cE, a, N)
    ax.plot(ref_data[i][:,0], conc_PEn/conc_PEn.sum(), c=clist[i])
ax.set_xlabel('m')
ax.set_ylabel(r'p(PE$_m$)')
ax.set_title(r'Lys (N=%.1f; a=%.3f)'%(N,a))

plt.tight_layout()
plt.show()

