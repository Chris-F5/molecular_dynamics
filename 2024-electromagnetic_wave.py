import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy

dimensions = 100
iterations = 800
dx = 1e-3
dt = 1e-12
wave = np.zeros((dimensions,dimensions))
wave_last = np.zeros((dimensions,dimensions))
wave[dimensions//2][dimensions//2] = 1
images = np.zeros((iterations+1, dimensions, dimensions))
c=2.998e8

kernel = [[0,  1, 0],
          [1, -4, 1],
          [0,  1, 0]]

images[0] = wave
for i in range(iterations):
  new_wave = 2 * wave - wave_last + ((c*dt/dx)**2) * scipy.ndimage.convolve(wave, kernel, mode='constant')
  #new_wave = np.maximum(new_wave, 0)
  wave_last = wave
  wave = new_wave
  images[i+1] = wave


fig,ax=plt.subplots()
ax.imshow(images[0],vmin=0.0, vmax=1.0)
artists = []
for i in range(0, len(images), 4):
  im = ax.imshow(images[i], vmin=0.0, vmax=0.5, animated=True)
  artists.append([im])
ani = animation.ArtistAnimation(fig, artists, interval=10, blit=True)
ani.save('./wave.gif',fps=24)

plt.show()
