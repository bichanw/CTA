# load mat
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression

# reading file content
mat_contents = sio.loadmat('mat/tmp.mat')
X = mat_contents['X']
Y = mat_contents['Y'].squeeze()
penalty = mat_contents['ops'][0]['mnr'][0]['penalty'][0][0][0]
C = 1 / float(mat_contents['ops'][0]['mnr'][0]['lambda'][0][0][0][0])

# fit model
try:
	clf = LogisticRegression(penalty=penalty,
							 C=C,
							 max_iter=1e4,
							 solver='liblinear').fit(X,Y)
except ValueError as ME:
	print('Error ')
sio.savemat('mat/pythonsave.mat',{'coef': clf.coef_.T, 'intercept': clf.intercept_})