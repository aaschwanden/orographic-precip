#!/usr/bin/env python
# Copyright (C) 2015 Andy Aschwanden
#

from argparse import ArgumentParser
import numpy as np
import pylab as plt

try:
    import pypismtools.pypismtools as ppt
except:
    import pypismtools as ppt

# Set up the option parser
parser = ArgumentParser()
parser.description = "Under construction."
parser.add_argument("-p", "--print_size", dest="print_mode",
                    choices=[
                        'onecol',
                        'medium',
                        'twocol',
                        'height',
                        'presentation',
                        'small_font',
                        'large_font',
                        '50mm',
                        '72mm'],
                    help="sets figure size and font size, available options are: \
                    'onecol','medium','twocol','presentation'", default="twocol")

options = parser.parse_args()
print_mode = options.print_mode

lw, pad_inches = ppt.set_mode(print_mode, aspect_ratio=0.6)

x1, zi = np.loadtxt('basal_topography.txt', unpack=True)
mask = np.zeros_like(zi)
mask[zi==-2e9] = 1
z = np.ma.array(data=zi, mask=mask)
x2, hi = np.loadtxt('surface_topography.txt', unpack=True)
mask = np.zeros_like(hi)
mask[hi==-2e9] = 1
h = np.ma.array(data=hi, mask=mask)

fig = plt.figure()
axU = fig.add_subplot(211)
axU.plot(x1 / 1e3, z, color='k', linewidth=2)
axU.set_ylabel('surface mass balance (m/yr)')

# Hide the right and top spines
axU.spines['right'].set_visible(False)
axU.spines['top'].set_visible(False)
axU.spines['bottom'].set_visible(False)
# Only show ticks on the left and bottom spines
axU.yaxis.set_ticks_position('left')
axU.tick_params(axis='x', bottom='off', top='off')
axU.set_xticks([])

axL = fig.add_subplot(212)
axL.plot(x1 / 1e3, z, color='k', linewidth=2)
axL.plot(x2 / 1e3, h, color='b', linewidth=2)
axL.set_xlabel('distance along profile (m)')
axL.set_ylabel('altitude (m a.s.l.)')

# Hide the right and top spines
axL.spines['right'].set_visible(False)
axL.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
axL.yaxis.set_ticks_position('left')
axL.tick_params(axis='x', top='off')
plt.tight_layout()
plt.savefig('ap_profile.pdf')
plt.close()
