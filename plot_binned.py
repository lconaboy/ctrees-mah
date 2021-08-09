import numpy as np
import matplotlib.pyplot as plt

z = np.loadtxt('./out_z.txt')
m = np.loadtxt('./out_mba.txt')
Mz = np.loadtxt('./out_mz.txt')
sz = np.loadtxt('./out_sz.txt')

cols = plt.cm.viridis(np.linspace(0., 1., num=len(m)))

fig, ax = plt.subplots(figsize=(6, 6))
for i in range(Mz.shape[0]):
    ax.errorbar(z, Mz[i, :], yerr=sz[i, :], label='{0:5.3f} h$^{{-1}}$ M$_\\odot$'.format(m[i]), color=cols[i])

ax.set_ylabel('M(z)/M(z={0:.2f})'.format(z[0]))
ax.set_xlabel('z')
ax.legend()
fig.savefig('out_binned.pdf', bbox_inches='tight')
