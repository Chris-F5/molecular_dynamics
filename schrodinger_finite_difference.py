import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy

iterations = 120000
dimensions = 100
save_interval = 200
dx = 3e-3
dt = 1e-7
#mass = 9.1e-31
mass = 1
#planck_constant = 6.626e-34
reduced_planck_constant = 1

wave_function = np.zeros(dimensions, dtype=np.complex128)
#wave_function[dimensions//2] = 1
wave_function[20:-20] = np.sin(np.linspace(0, np.pi, len(wave_function[20:-20])))
wave_function[20:-20] += np.sin(np.linspace(0, 2*np.pi, len(wave_function[20:-20])))
wave_function = wave_function / np.sum(np.absolute(wave_function)**2)*dx
#potential = (5e5)*(0-np.sin(np.linspace(0, np.pi, dimensions)))*0

# The potentials are behaving strangely. Seems like the wave likes to stay in
# a high potential well.
potential = np.zeros(dimensions)
potential[20:-20] = 1e2 # 1e6
#potential[:20] = 1e6
#potential[-20:] = 1e6
#potential = np.linspace(0, 1e10, dimensions)
#wave_function = np.sin(np.linspace(0, np.pi, dimensions))
kernel = [1, -2, 1]

image = np.zeros((iterations//save_interval, dimensions), dtype=np.complex128)
for i in range(iterations):
  if i % save_interval == 0:
    image[i//save_interval] = wave_function
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
  wave_function = wave_function / np.sum(np.absolute(wave_function)**2)*dx

fig,(ax,ax2) = plt.subplots(2)
ax.imshow(np.absolute(image.T)**2)
ax2.plot(np.linspace(0, dx*dimensions, dimensions), np.absolute(image[0])**2)

artists = []
for i in range(0, iterations//save_interval):
  l = ax2.plot(np.linspace(0, dx*dimensions, dimensions), np.absolute(image[i])**2, c='black')
  artists.append(l)

ani = animation.ArtistAnimation(fig, artists, interval=50, blit=True)

plt.show()
