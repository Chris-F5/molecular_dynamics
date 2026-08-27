import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy
import math

duration = 4
hbar = 1
m = 1
a = 6
dx = 0.1
dt = 0.0001
k = np.pi/a
t = 0
c1 = math.sqrt(3/10)
c2 = math.sqrt(7/10)

dimensions = int(a/dx)

# Amazingly, the two following lines which each set x in a slightly different
# way result in tremendous differences in simulation outcome. The first of the
# two produces a stable simulation that behaves as expected. The second results
# in a simulation which becomes unstable quickly. Clearly, such a sensitive
# system is not desirable. I should find a way to make it more resiliant.
#x = np.linspace(0, a, dimensions + 2)[1:-1]
x = np.linspace(0+dx, a-dx, dimensions)


print(x)

wave_function = (c1*np.sin(x*k) + c2*np.sin(2*x*k)).astype(np.complex128)
V = np.zeros(dimensions)
PP = np.square(np.absolute(wave_function))

anim_dt = 5e-2
anim_index = 0
anim_frames = int(duration/anim_dt)
anim_image = np.zeros((anim_frames, dimensions), dtype=np.complex128)

while t<duration:
  if t >= anim_index * anim_dt:
    #anim_image[anim_index] = PP
    anim_image[anim_index] = wave_function
    anim_index += 1
    #print(f"{anim_index}/{anim_frames}")

  wave_function \
    = wave_function \
    + (1j * hbar * dt / (2*m*(dx**2))) * scipy.ndimage.convolve(wave_function, [1,-2,1], mode='constant') \
    - 1j * (dt/hbar) * V * wave_function
  
  wave_function = wave_function / math.sqrt(np.sum(np.square(np.absolute(wave_function)))*dx)

  t = t + dt

fig,ax = plt.subplots()
artists = []
for i in range(anim_frames):
  probability = ax.plot(x, np.square(np.absolute(anim_image[i])), color='black')
  real = ax.plot(x, np.real(anim_image[i]), color='red')
  imaginary = ax.plot(x, np.imag(anim_image[i]), color='blue')
  artists.append(probability + real + imaginary)
anim = animation.ArtistAnimation(fig, artists, interval=100)
plt.show()
