from PIL import Image
import csv
import numpy as np

img = Image.open("test.bmp")
pixels = img.load()
data = np.zeros((img.size[0],img.size[1]))
for i in range(img.size[0]):    # for every pixel:
    for j in range(img.size[1]):
        data[i,j] = int(pixels[i, j])
np.savetxt("cleanData.csv", data, fmt='%i', delimiter=",")
print(img.size[0])
print(img.size[1])

