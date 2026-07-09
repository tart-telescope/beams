"""Load measured/simulated beam samples and fit a :class:`Beam`.

Input format is a JSON list of records::

    [{"el": 90.0, "az": 0.0, "gain": 1.0},
     {"el": 45.0, "az": 30.0, "gain": 0.5}, ...]

``el`` is elevation above the horizon and ``az`` is azimuth, both in **degrees**.
The beam is assumed to point at the **zenith**, so a sample at elevation ``el``
and azimuth ``az`` corresponds to the sky unit vector

    x = cos(el) cos(az),  y = cos(el) sin(az),  z = sin(el)

In the zenith-pointing beam frame this gives ``w = z = sin(el)``: ``el = 90`` is
boresight (``w = 1``), ``el = 0`` is the horizon (``w = 0``), and ``el < 0`` is
behind the antenna (``w < 0``, dropped by the fit).
"""

# Copyright (c) 2026 Tim Molteno (tim@elec.ac.nz)

import json

import numpy as np

from .beam import Beam

ZENITH = (0.0, 0.0, 1.0)


def elaz_to_vec(el_deg, az_deg):
    """Convert elevation/azimuth (degrees) to unit sky vectors, shape ``(...,3)``.

    Assumes a zenith-pointing frame: ``z`` is up, ``az`` is measured from the
    ``x`` axis toward ``y``.
    """
    el = np.radians(np.asarray(el_deg, dtype=float))
    az = np.radians(np.asarray(az_deg, dtype=float))
    cos_el = np.cos(el)
    return np.stack([cos_el * np.cos(az), cos_el * np.sin(az), np.sin(el)],
                    axis=-1)


def samples_from_records(records):
    """Turn a list of ``{el, az, gain}`` dicts into ``(s_hat, gain)`` arrays."""
    if not records:
        raise ValueError("no records to load")
    el = np.array([r["el"] for r in records], dtype=float)
    az = np.array([r["az"] for r in records], dtype=float)
    gain = np.array([r["gain"] for r in records], dtype=float)
    return elaz_to_vec(el, az), gain


def load_elaz_json(path):
    """Read an el/az/gain JSON file and return ``(s_hat, gain)`` arrays."""
    with open(path) as f:
        records = json.load(f)
    return samples_from_records(records)


def fit_from_json(path, degree=8, q=2, ridge=1e-6, weights=None, **kw):
    """Load an el/az/gain JSON file and fit a zenith-pointing :class:`Beam`.

    Extra keyword arguments (``name``, ``frequency_hz``, ...) are passed to
    :meth:`Beam.fit`. The beam can be re-pointed afterwards with
    :meth:`Beam.set_pointing`.
    """
    s_hat, gain = load_elaz_json(path)
    return Beam.fit(s_hat, gain, degree=degree, q=q, ridge=ridge,
                    weights=weights, boresight=ZENITH, **kw)
