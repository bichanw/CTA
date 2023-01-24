# load mat
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression

# reading
# change this to input
mat_contents = sio.loadmat('mat/tmp.mat')
X = mat_contents['X']
Y = mat_contents['Y'].squeeze()

# fit model
clf = LogisticRegression().fit(X,Y)
sio.savemat('mat/pythonsave.mat',{'coef': clf.coef_.T, 'intercept': clf.intercept_})