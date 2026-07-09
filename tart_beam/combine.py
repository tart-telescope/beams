"""Combine overlapping beams into full-sphere HEALPix products.

See ``DESIGN.md`` section 7. ``healpy`` is imported lazily so the core beam
model has no hard dependency on it.
"""

# Copyright (c) 2026 Tim Molteno (tim@elec.ac.nz)

import numpy as np


def healpix_directions(nside, nest=False):
    """Unit direction vectors for every HEALPix pixel, shape ``(npix, 3)``."""
    import healpy as hp

    npix = hp.nside2npix(nside)
    x, y, z = hp.pix2vec(nside, np.arange(npix), nest=nest)
    return np.stack([x, y, z], axis=-1)


def coverage_map(beams, nside, nest=False):
    """Summed power response ``A = sum_k |B_k|**2`` over the full sphere."""
    s_hat = healpix_directions(nside, nest)
    total = np.zeros(s_hat.shape[0])
    for beam in beams:
        total += beam.evaluate(s_hat) ** 2
    return total


def partition_of_unity(beams, nside, nest=False, eps=0.0):
    """Per-beam blending weights ``w_k = B_k**2 / sum_j B_j**2``.

    Returns an array of shape ``(n_beams, npix)``. Pixels with no coverage
    (denominator ``<= eps``) get weight zero across all beams.
    """
    s_hat = healpix_directions(nside, nest)
    powers = np.stack([beam.evaluate(s_hat) ** 2 for beam in beams], axis=0)
    denom = powers.sum(axis=0)
    covered = denom > eps
    weights = np.zeros_like(powers)
    weights[:, covered] = powers[:, covered] / denom[covered]
    return weights


def mosaic(beams, images, sigmas=None, nside=None, nest=False):
    """Optimal inverse-variance mosaic of per-beam sky estimates.

        I = sum_k B_k I_k / sigma_k**2  /  sum_k B_k**2 / sigma_k**2

    Parameters
    ----------
    beams : list[Beam]
    images : list[array]
        Per-beam HEALPix maps ``I_k`` (same nside as the beams are sampled on).
    sigmas : list[float] or None
        Per-beam noise sigma; defaults to all ones.
    nside : int
        HEALPix resolution; inferred from ``images`` if omitted.

    Returns the combined map plus the (effective-weight) denominator, so callers
    can mask uncovered pixels.
    """
    import healpy as hp

    images = [np.asarray(im, dtype=float) for im in images]
    if nside is None:
        nside = hp.npix2nside(images[0].size)
    if sigmas is None:
        sigmas = np.ones(len(beams))

    s_hat = healpix_directions(nside, nest)
    num = np.zeros(images[0].shape)
    den = np.zeros(images[0].shape)
    for beam, im, sig in zip(beams, images, sigmas):
        b = beam.evaluate(s_hat)
        inv_var = 1.0 / (sig * sig)
        num += b * im * inv_var
        den += b * b * inv_var

    out = np.full_like(num, hp.UNSEEN)
    good = den > 0.0
    out[good] = num[good] / den[good]
    return out, den
