import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as mpl
import scipy

size = 1
dimensions = 100
duration = 1e-1
dx = size / dimensions
dt = 1e-3
iterations = int(duration/dt)
concentration = np.zeros((dimensions, dimensions))
concentration[dimensions//10][dimensions//10] = 1
kernel = [[0,  1, 0],
          [1, -4, 1],
          [0,  1, 0]]
diffusion_constant = 1e-3
image = np.zeros((iterations + 1, dimensions, dimensions))
image[0] = concentration

for i in range(iterations):
  laplacian = (2.0 / (dx**2)) * scipy.ndimage.convolve(concentration, kernel, mode='wrap')
  concentration = concentration + diffusion_constant *  dt * laplacian
  image[i+1] = concentration

fig,ax = plt.subplots()
ax.imshow(image[0], vmin=0.0, vmax=1.0)
artists = []
for i in range(0, len(image)):
  im = ax.imshow(image[i], vmin=0.0, vmax=1.0, animated=True)
  txt = ax.text(10, 10, f"{i}/{len(image)}", c='green')
  artists.append([im, txt])
print(artists[-1])
ani = animation.ArtistAnimation(fig, artists, interval=1, blit=True, repeat_delay=1000)

plt.show()
