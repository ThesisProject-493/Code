from PIL import Image
import csv
import numpy as np

data = np.loadtxt("restoredData40_4.csv",delimiter=",")
data = np.transpose(data)
result = Image.fromarray(data.astype(np.uint8))
result.save('restoredImage40_4.bmp')