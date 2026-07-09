"""End-to-end tests on a real HEALPix sphere (skipped if healpy is absent)."""

# Copyright (c) 2026 Tim Molteno (tim@elec.ac.nz)

import numpy as np
import pytest

hp = pytest.importorskip("healpy")

from tart_beam import Beam, coverage_map, healpix_directions, mosaic, partition_of_unity


@pytest.fixture
def tiled_beams():
    """A few overlapping zenith-ish beams pointed around the upper sky."""
    idx, coef = [(0, 0), (2, 0)], [1.0, -0.3]

    def unit(theta_deg, phi_deg):
        t, p = np.radians(theta_deg), np.radians(phi_deg)
        return [np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)]

    pointings = [unit(0, 0), unit(50, 0), unit(50, 120), unit(50, 240)]
    return [Beam(idx, coef, boresight=p, q=2) for p in pointings]


def test_healpix_directions_are_unit():
    v = healpix_directions(nside=8)
    assert v.shape == (hp.nside2npix(8), 3)
    assert np.allclose(np.linalg.norm(v, axis=-1), 1.0)


def test_coverage_nonnegative_and_partial(tiled_beams):
    cov = coverage_map(tiled_beams, nside=16)
    assert cov.shape == (hp.nside2npix(16),)
    assert np.all(cov >= 0.0)
    # upper hemisphere is covered; the deep south is not -> partial sky
    frac = (cov > 0).mean()
    assert 0.0 < frac < 1.0


def test_partition_of_unity_on_healpix(tiled_beams):
    nside = 16
    w = partition_of_unity(tiled_beams, nside=nside)
    assert w.shape == (len(tiled_beams), hp.nside2npix(nside))
    colsum = w.sum(axis=0)
    cov = coverage_map(tiled_beams, nside=nside)
    covered = cov > 0
    assert np.allclose(colsum[covered], 1.0)
    assert np.allclose(colsum[~covered], 0.0)


def test_mosaic_recovers_sky_on_healpix(tiled_beams):
    nside = 16
    s_hat = healpix_directions(nside)
    sky = 3.5
    # apparent per-beam images = beam * true sky
    imgs = [b.evaluate(s_hat) * sky for b in tiled_beams]
    out, den = mosaic(tiled_beams, imgs, nside=nside)

    covered = den > 0
    assert np.allclose(out[covered], sky, atol=1e-6)
    # uncovered pixels are flagged UNSEEN
    assert np.all(out[~covered] == hp.UNSEEN)
