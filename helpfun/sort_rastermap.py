import numpy as np
import matplotlib.pyplot as plt
# importing rastermap
# (this will be slow the first time since it is compiling the numba functions)
from rastermap import Rastermap, utils
from scipy.stats import zscore

# download processed hippocampus recording from Grosmark & Buzsaki 2016
spks = np.load('/mnt/cup/people/bichanw/SpikeSorting/Codes/CTA/mat/spks.npy')

# spks is neurons by time
# (each timepoint is 200 ms)
# spks = dat["spks"]
n_neurons, n_time = spks.shape
print(f"{n_neurons} neurons by {n_time} timepoints")
# zscore activity (each neuron activity trace is then mean 0 and standard-deviation 1)
spks = zscore(spks, axis=1)

# location of the rat and speed
# loc2d = dat["loc2d"] # 2D location
# loc_signed = dat["loc_signed"] # left runs are positive and right runs are negative
# speed = (np.diff(loc2d, axis=0)**2).sum(axis=1)**0.5
# speed = np.concatenate((np.zeros((1,)), speed), axis=0)

# # which neurons in the recording are pyramidal cells
# pyr_cells = dat["pyr_cells"].astype("int")

model = Rastermap(n_clusters=None, # None turns off clustering and sorts single neurons 
                  n_PCs=64, # use fewer PCs than neurons
                  locality=0.15, # some locality in sorting (this is a value from 0-1)
                  time_lag_window=15, # use future timepoints to compute correlation
                  grid_upsample=0, # 0 turns off upsampling since we're using single neurons
                ).fit(spks)
y = model.embedding # neurons x 1
isort = model.isort
np.save('mat/sorted_rastermap.npy', isort)

# ax = plt.subplot()
# # xmin = 0
# # xmax = xmin + 2000
# ax.imshow(spks[isort, :2000], cmap="gray_r", vmin=0, vmax=1.2, aspect="auto")
# ax.set_xlabel("time")
# ax.set_ylabel("superneurons")
# plt.savefig('tmp.png', dpi=300)



