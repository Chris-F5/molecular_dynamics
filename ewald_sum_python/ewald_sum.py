import math
import numpy as np
import scipy
from hypothesis import given, example, strategies as st
import matplotlib.pyplot as plt


def generate_periodic_images(image_radius):
    r = image_radius
    x, y, z = map(np.ravel, np.meshgrid(*[np.arange(-r, r+1) for i in range(3)]))
    images = np.column_stack([x, y, z])
    return images


def generate_pair_list(particle_count, image_radius):
    images = generate_periodic_images(image_radius)
    i, j, image_index = map(np.ravel, np.meshgrid(
        np.arange(particle_count), np.arange(particle_count),
        np.arange(len(images))
        ))
    n = images[image_index]
    self_interaction = (i == j) & np.all(n == [0, 0, 0], axis=1)
    i, j, n = (a[~self_interaction] for a in (i, j, n))
    return i, j, n


@given(r=st.integers(0, 5))
def test_generate_periodic_images(r):
    images = generate_periodic_images(r)
    assert len(images) == (r*2+1)**3, "invalid count"
    assert len(np.unique(images, axis=0)) == len(images), "contains duplicates"
    assert [0, 0, 0] in images.tolist(), "does not contain zero"


@given(N=st.integers(0, 10), r=st.integers(0, 3))
@example(N=0, r=2)
@example(N=2, r=0)
def test_generate_pair_list(N, r):
    i, j, n = generate_pair_list(N, r)
    assert i.ndim == 1, "invalid shape"
    assert j.ndim == 1, "invalid shape"
    assert n.ndim == 2, "invalid shape"
    assert n.shape[1] == 3, "invalid shape"
    assert len(i) == len(j) == len(n), "inconsistent sizes"
    assert len(i) == N*N*((r*2+1)**3)-N, "incorrect count"
    assert np.none(i == 0 & j == 0 & np.all(n == [0, 0, 0], axis=1))


if __name__ == "__main__":
    np.random.seed(123)
    N = 2
    L = 1
    assert N % 2 == 0
    position = np.random.random((N, 3)) * L
    charge = np.tile([1, -1], N//2)

    def smear_potential(image_radius, alpha):
        k = (2*np.pi/L) * generate_periodic_images(image_radius)
        k = k[~np.all(k == [0, 0, 0], axis=1)]
        kdotr = np.einsum("ik,jk->ij", k, position)
        rhok = np.einsum("i,ki->k", charge, np.exp(1j*kdotr))
        rhoksquared = np.square(np.absolute(rhok))
        ksquared = np.square(np.linalg.norm(k, axis=1))
        potential = 2*np.pi / (L**3) * np.sum(
            rhoksquared * np.exp(-ksquared/(4*alpha)) / ksquared
        )
        return potential

    def self_interaction_potential(alpha):
        return -np.sqrt(alpha / np.pi) * np.sum(np.square(charge))

    def real_space_potential(alpha, cutoff_radius):
        i, j, n = generate_pair_list(N, math.ceil(cutoff_radius/L))
        offset = position[i] - (position[j] + n * L)
        distance = np.linalg.norm(offset, axis=1)
        in_cutoff = distance < cutoff_radius
        i, j, n, offset, distance = (a[in_cutoff] for a in (
            i, j, n, offset, distance))
        potential = 1/2 * np.sum(
            charge[i] * charge[j] * scipy.special.erfc(np.sqrt(alpha)*distance) / distance)
        return potential

    #fig, axs = plt.subplots(1, 2, sharey=True)
    #for alpha in np.linspace(5, 50, 10):
    #    k_radius = np.arange(1, 20)
    #    smear_p = np.array([smear_potential(radius, alpha) for radius in k_radius])
    #    r_radius = np.linspace(0.1, 5)
    #    real_p = np.array([real_space_potential(alpha, radius) for radius in r_radius])
    #    axs[0].plot(k_radius, smear_p, label=f"a={alpha},smear")
    #    axs[1].plot(r_radius, real_p, label=f"a={alpha},real")
    for alpha in np.linspace(1, 5, 10):
        radius = np.arange(1, 30)
        smear_p = np.array([smear_potential(r, alpha) for r in radius])
        real_p = np.array([real_space_potential(alpha, r/4) for r in radius])
        self_p = self_interaction_potential(alpha)
        plt.plot(radius, smear_p + real_p + self_p, label=f"a={alpha}")
    plt.legend()
    plt.show()


    #cutoff_radius = 50
    #i, j, n = generate_pair_list(N, math.ceil(cutoff_radius/L))
    #pair_offset = position[i] - (position[j] + n * L)
    #pair_distance = np.linalg.norm(pair_offset, axis=1)
    #in_cutoff = pair_distance < cutoff_radius
    #i, j, n, pair_offset, pair_distance = (a[in_cutoff] for a in (
    #    i, j, n, pair_offset, pair_distance))
    #pair_potential = charge[i] * charge[j] / pair_distance

    #pair_sort = np.argsort(pair_distance)
    #distance = pair_distance[pair_sort]
    #potential = np.cumsum(pair_potential[pair_sort])

    #plt.plot(distance, potential)
    #plt.legend()
    #plt.show()
