"""Tests for full-sphere combination.

healpy is not required: we stub ``healpix_directions`` with a small fixed set of
sky vectors, and inject a minimal fake ``healpy`` module for ``mosaic``.
"""

import sys
import types

import numpy as np
import pytest

import tart_beam.combine as combine
from tart_beam import Beam

UNSEEN = -1.6375e30


@pytest.fixture
def fixed_pixels(monkeypatch):
    """Four directions: +z, -z (gap behind both), +x, and a tilt."""
    pts = np.array([
        [0.0, 0.0, 1.0],    # zenith   -> only the z-beam
        [0.0, 0.0, -1.0],   # nadir    -> behind both -> gap
        [1.0, 0.0, 0.0],    # horizon  -> only the x-beam
        [0.6, 0.0, 0.8],    # tilted   -> both beams see it
    ])
    pts = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    monkeypatch.setattr(combine, "healpix_directions",
                        lambda nside, nest=False: pts)
    return pts


@pytest.fixture
def two_beams():
    idx, coef = [(0, 0)], [1.0]
    return [Beam(idx, coef, boresight=[0, 0, 1], q=2),
            Beam(idx, coef, boresight=[1, 0, 0], q=2)]


def test_coverage_map(fixed_pixels, two_beams):
    cov = combine.coverage_map(two_beams, nside=1)
    y00 = 1.0 / np.sqrt(4 * np.pi)
    # zenith: only beam-z at full response (w=1, taper=1) -> b = y00
    assert np.isclose(cov[0], y00 ** 2)
    # nadir: behind both antennas -> zero
    assert cov[1] == 0.0
    # tilted pixel: both beams see it, so it equals the sum of their powers
    bz, bx = (b.evaluate(fixed_pixels[3][None])[0] for b in two_beams)
    assert np.isclose(cov[3], bz ** 2 + bx ** 2)
    assert cov[3] > 0.0


def test_partition_of_unity_sums_to_one(fixed_pixels, two_beams):
    w = combine.partition_of_unity(two_beams, nside=1)
    assert w.shape == (2, 4)
    colsum = w.sum(axis=0)
    covered = np.array([True, False, True, True])
    assert np.allclose(colsum[covered], 1.0)
    assert colsum[1] == 0.0  # the gap gets zero weight everywhere
    assert np.all(w >= 0.0)


def _fake_healpy():
    return types.SimpleNamespace(UNSEEN=UNSEEN, npix2nside=lambda n: 1)


def test_mosaic_recovers_true_sky(fixed_pixels, two_beams, monkeypatch):
    """Apparent images I_k = B_k * sky -> the mosaic recovers `sky`."""
    monkeypatch.setitem(sys.modules, "healpy", _fake_healpy())
    sky = 7.0
    bz = two_beams[0].evaluate(fixed_pixels)
    bx = two_beams[1].evaluate(fixed_pixels)
    imgs = [bz * sky, bx * sky]  # beam-attenuated (apparent) per-beam images
    out, den = combine.mosaic(two_beams, imgs, sigmas=[1.0, 1.0], nside=1)

    # wherever there is any coverage, the true sky value is recovered
    covered = np.array([True, False, True, True])
    assert np.allclose(out[covered], sky)
    # nadir gap is masked
    assert out[1] == UNSEEN
    assert den[1] == 0.0


def test_mosaic_noise_weighting(fixed_pixels, two_beams, monkeypatch):
    """A noisier beam is down-weighted toward the cleaner beam's estimate."""
    monkeypatch.setitem(sys.modules, "healpy", _fake_healpy())
    bz = two_beams[0].evaluate(fixed_pixels)
    bx = two_beams[1].evaluate(fixed_pixels)
    # apparent images implying different true skies: 4 (z-beam) vs 10 (x-beam)
    imgs = [bz * 4.0, bx * 10.0]
    out_eq, _ = combine.mosaic(two_beams, imgs, sigmas=[1.0, 1.0], nside=1)
    out_znoisy, _ = combine.mosaic(two_beams, imgs, sigmas=[10.0, 1.0], nside=1)
    # at the tilted (overlap) pixel, distrusting the z-beam pulls toward 10
    assert out_znoisy[3] > out_eq[3]
