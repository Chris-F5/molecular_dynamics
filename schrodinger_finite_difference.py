import numpy as np
import matplotlib.pyplot as plt
import scipy

iterations = 30000
dimensions = 100
dx = 3e-3
dt = 1e-7
#mass = 9.1e-31
mass = 1
#planck_constant = 6.626e-34
reduced_planck_constant = 1

wave_function = np.zeros(dimensions, dtype=np.complex128)
potential = (5e5)*(0-np.sin(np.linspace(0, np.pi, dimensions)))*0
potential[:dimensions//3] = 1e6
print(potential)
#potential = np.linspace(0, 1e10, dimensions)
wave_function[dimensions//2] = 1
#wave_function = np.sin(np.linspace(0, np.pi, dimensions))
wave_function = wave_function / np.sum(np.absolute(wave_function))
kernel = [1, -2, 1]

image = np.zeros((iterations+1, dimensions), dtype=np.complex128)
image[0] = wave_function
for i in range(iterations):
  wave_function \
    = wave_function \
    + (1j * reduced_planck_constant * dt / (2*mass*(dx**2))) * scipy.ndimage.convolve(wave_function, kernel, mode='constant') \
    - 1j * (dt/reduced_planck_constant) * potential * wave_function
  # The evolution function we use massivly increases the sum of the absolutes,
  # so we must make the following normalisation to ensure the square of the
  # wave function can still be interpreted as probability density. However, this
  # causes problems if a particular cell gets some massive number in it (like
  # what might happen at a sharp potential boundry). This makes if difficult for
  # this numerical scheme to handle scenarios with steep potential boundries.
  wave_function = wave_function / np.sum(np.absolute(wave_function))
  image[i+1] = wave_function

height = 200
reduced_image = np.zeros((height, dimensions))
for i in range(height):
  reduced_image[i] = image[i*(iterations//height)]

print(image)
fig,ax = plt.subplots()
ax.imshow(np.absolute(reduced_image), vmax=0.1)
plt.show()
