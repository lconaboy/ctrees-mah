"""Really nasty script to plot the output of mah.c.

TODO
   - errorbars

"""

import numpy as np
import matplotlib.pyplot as plt
from plot import get_cols


def minmax(x):
    """Return min and max of array x, allowing for NaN values"""
    if np.any(np.isnan(x)):
        mask = np.logical_not(np.isnan(x))
    else:
        mask = np.ones_like(x)
        
    return np.array([np.min(x[mask]), np.max(x[mask])])


def nan_minmax(x):
    """Turns out numpy already has this"""
    return np.array([np.nanmin(x), np.nanmax(x)])


def get_bins(x):
    """Returns a set of bins going up decades in mass as 1.0En ->
    1.0En+1, smallest and largest n is determined by the min and max
    final values of x.

    :param x: (array) data to bin

    :returns: bin edges

    :rtype: (array)
    """
    
    mi, ma = nan_minmax(x[:, 0])  # bin based on final mass
    bin_edges = 10.**np.arange(int(np.log10(mi)), int(np.log10(ma)) + 2)

    return bin_edges


def bin_mah(x):
    """Bins the mass accretion histories, x, in the format returned by
    read_mah.  Calculates the binned MAH in whole decades (i.e. 10^5
    -> 10^6 Msol/h) binned by the final halo mass and the normalised
    (i.e. M(z)/M(z_f)) MAH binned in the same way.

    :param x: (array) masses in format returned by read_mah

    :returns: binned MAH, error on binned MAH, binned normalised MAH,
              error on binned normalised MAH, bin edges

    :rtype: (array, array, array, array, array)
    """
    bin_edges = get_bins(x)
    nbins = len(bin_edges) - 1

    # Preallocate arrays for mean and standard deviation
    avg = np.zeros((nbins, x.shape[1]))
    std = np.zeros((nbins, x.shape[1]))

    norm_avg = np.zeros((nbins, x.shape[1]))
    norm_std = np.zeros((nbins, x.shape[1]))

    for i in range(nbins):
        # Get a mask containing the indices of the data in each bin
        mask = np.logical_and(x[:, 0] >= bin_edges[i],
                              x[:, 0] < bin_edges[i+1])

        # Get the actual data in the mask
        tmp = x[mask, :]
        avg[i, :] = np.nanmean(tmp, axis=0)
        std[i, :] = np.nanstd(tmp, axis=0)

        # Normlise to final mass
        tmp = norm_mvir(tmp)
        norm_avg[i, :] = np.nanmean(tmp, axis=0)
        norm_std[i, :] = np.nanstd(tmp, axis=0)

        
        # Check distribution of mvir vals
        # for i in range(x.shape[1]):
        #     mi, ma = nan_minmax(x[mask, i])
        #     bin_edges = 10.**np.linspace(np.log10(mi), np.log10(ma),
        #                                  num=51)
        #     plt.figure()
        #     plt.hist(x[mask, i], bins=bin_edges)
        #     plt.xscale('log')
        #     plt.show()
    
    return avg, std, norm_avg, norm_std, bin_edges


def bin_norm_mah(x):
    return


def norm_mvir(x):
    """Normalises the masses to be a function of the final snapshot
    mass"""

    # python modifies mutable objects, so need to copy
    y = np.zeros_like(x)
    
    for i in range(x.shape[0]):
        y[i, :] = x[i, :] / x[i, 0]
        print('final halo mass: {0:.2e} Msol/h'.format(x[i, 0]))
    return y


def read_mah(fn='out.txt'):
    """Reads in the out.txt files produced by mah.c. Expects the file
    format:

    # nroots nsnaps #
    ...
    #
    ...

    where '...' corresponds to the scale factor and progenitor mass,
    one scale factor per line, later time first.

    """

    with open(fn, 'r') as f:
        i = 0
        j = 0

        # Grab number of trees and snaps
        l = f.readline().strip('\n').split()
        d = (int(l[1]), int(l[2]))

        aexp = np.zeros(d)
        mvir = np.zeros(d)
        # z = np.zeros(d)
        aexp[:] = np.nan
        mvir[:] = np.nan
        # z[:] = np.nan

        for l in f:
            if '#' in l:
                # if i % 10 == 0: ax.plot(z[i, :], mvir[i, :], c='k', alpha=0.5)
                i += 1
                j = 0
            else:
                l = l.strip('\n').split(' ')
                aexp[i, j] = float(l[0])
                mvir[i, j] = float(l[1])
                # z[i, j] = 1./aexp[i, j] - 1.  # slow, could do this outside loop
                j += 1

    z = (1. / aexp) - 1.
    # Some scale factors would be NaN, so we can get all of the
    # redshifts by using unique on z -- unique also handily sorts the
    # array.
    all_z = np.unique(z)
    
    return z, all_z, mvir


z, all_z, mvir = read_mah()
# norm_mvir = norm_mvir(mvir)

# fig, ax = plt.subplots(figsize=(6, 6))
# ax.loglog(all_z, mvir)
# ax.set_xlabel('z')
# ax.set_ylabel('M (M$_\odot$/h)')
# ax.set_yscale('log')
# fig.savefig('./out.pdf')

avg, std, norm_avg, norm_std, bin_edges = bin_mah(mvir)
nbins = avg.shape[0]
nz = avg.shape[1]
c = get_cols(nbins)

fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlabel('z')
# ax.set_ylabel('M (M$_\odot$/h)')
ax.set_ylabel('M(z) (M$_\odot$/h)')
ax.set_yscale('log')
ax.set_xscale('log')
for i in range(nbins):
    l = '{0:.2e} $\leq$ M (M$_\odot$/h) $<$ {1:.2e}'.format(bin_edges[i],
                                                            bin_edges[i+1])
    ax.plot(all_z[0:nz], avg[i, :], c=c[i], label=l)
    # ax.errorbar(all_z[0:avg.shape[1]], avg[i, :], yerr=std[i, :],
    #             c=c[i], label=l)
ax.legend()
fig.tight_layout()
fig.savefig('out_binned.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlabel('z')
# ax.set_ylabel('M (M$_\odot$/h)')
ax.set_ylabel('M(z)/M({0:.2f})'.format(all_z[0]))
ax.set_yscale('log')
ax.set_xscale('log')
for i in range(nbins):
    l = '{0:.2e} $\leq$ M (M$_\odot$/h) $<$ {1:.2e}'.format(bin_edges[i],
                                                            bin_edges[i+1])
    ax.plot(all_z[0:nz], norm_avg[i, :], c=c[i], label=l)
    # ax.errorbar(all_z[0:avg.shape[1]], norm_avg[i, :], yerr=norm_std[i, :],
    #             c=c[i], label=l)
ax.legend()
fig.tight_layout()
fig.savefig('out_norm_binned.pdf', bbox_inches='tight')
