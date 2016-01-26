from PIL import Image
import csv
import numpy as np

data = np.loadtxt("restoredData.csv",delimiter=",")
data = np.transpose(data)
result = Image.fromarray(data.astype(np.uint8))
result.save('restoredImage.bmp')