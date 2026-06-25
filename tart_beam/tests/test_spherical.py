"""Tests for the real spherical-harmonic basis."""

import numpy as np
import pytest

from tart_beam.spherical import basis_matrix, real_sph_harm, sh_indices


def test_sh_indices_count():
    # number of harmonics up to degree N is (N+1)**2
    for N in range(6):
        assert len(sh_indices(N)) == (N + 1) ** 2


def test_y00_is_constant():
    theta = np.array([0.0, 0.5, 1.0, 2.0])
    phi = np.array([0.0, 1.0, 2.0, 3.0])
    y00 = real_sph_harm(0, 0, theta, phi)
    assert np.allclose(y00, 1.0 / np.sqrt(4 * np.pi))


def test_orthonormality_on_grid():
    # numerically integrate Y_a * Y_b over the sphere; expect the identity
    indices = sh_indices(3)
    nt, npp = 200, 400
    theta = (np.arange(nt) + 0.5) * np.pi / nt
    phi = (np.arange(npp) + 0.5) * 2 * np.pi / npp
    T, P = np.meshgrid(theta, phi, indexing="ij")
    w = np.sin(T) * (np.pi / nt) * (2 * np.pi / npp)  # solid-angle weights

    M = basis_matrix(indices, T.ravel(), P.ravel())  # (npts, nbasis)
    gram = (M * w.ravel()[:, None]).T @ M
    assert np.allclose(gram, np.eye(len(indices)), atol=2e-3)
