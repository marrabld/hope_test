#!/usr/bin/env python

import numpy as np
import os, sys, string, re
import matplotlib.pyplot as plt


Common_Wavelength = np.array([410.0, 430.0, 450.0, 470.0, 490.0, 510.0, 530.0, 550.0, 570.0, 590.0, 610.0, 630.0, 650.0, 670.0, 690.0, 710.0, 730.], order='C', dtype=np.float)

HOPE_rrs_array  = np.load("HOPE_rrs_17Bands.npy")
ModelInputs_arr = np.load("ModelInputs_17Bands.npy")

Number_Spectra, Number_WaveL = HOPE_rrs_array.shape


for i in xrange(0, Number_Spectra):

   P = ModelInputs_arr[i, 0]
   G = ModelInputs_arr[i, 1]
   X = ModelInputs_arr[i, 2]
   H = ModelInputs_arr[i, 3]
   Y = ModelInputs_arr[i, 4]
   S = ModelInputs_arr[i, 5]

   Title_Str = "P = "+str(P)+"; G = "+str(G)+"; X = "+str(X)+"; H = "+str(H)+"; Y = "+str(Y)+"; S = "+str(S)

   plt.title(Title_Str)
   plt.plot(Common_Wavelength, HOPE_rrs_array[i])
plt.show()
