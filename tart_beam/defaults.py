"""Default / reference beam patterns for TART antennas.

The main entry point is :func:`base_tart_beam`, which returns the nominal primary
beam pattern of a TART antenna element: uniform above 10° elevation, tapering
smoothly to zero at the horizon.
"""

# Copyright (c) 2026 Tim Molteno (tim@elec.ac.nz)

import functools

import numpy as np

from .beam import Beam

__all__ = ["base_tart_beam"]

_ELEVATION_TRANSITION_DEG = 10.0
"""Elevation (degrees) above which the base beam is uniform."""


def _fib_front_hemisphere(n_total):
    """Fibonacci-sphere samples covering the front (w >= 0) hemisphere.

    Returns ``(s_hat, w, el_deg)``, each shaped ``(n_front,)`` or ``(n_front, 3)``.
    """
    i = np.arange(n_total) + 0.5
    z = 1.0 - 2.0 * i / n_total
    front = z >= 0
    z_f = z[front]
    n_f = np.count_nonzero(front)
    r = np.sqrt(np.clip(1.0 - z_f * z_f, 0.0, 1.0))
    phi = np.pi * (1.0 + np.sqrt(5.0)) * i[front]
    s_hat = np.stack([r * np.cos(phi), r * np.sin(phi), z_f], axis=-1)
    w = z_f
    el_deg = np.degrees(np.arcsin(np.clip(w, -1.0, 1.0)))
    return s_hat, w, el_deg


def _cosine_taper_profile(el_deg, transition_deg):
    """Smooth cosine taper from 1 at *transition_deg* to 0 at the horizon.

    For elevations >= *transition_deg* the return value is 1.
    The taper is applied in ``w = sin(el)`` space so the transition is linear in
    the direction cosine ``w`` — natural for the beam model's ``w**q`` factor.
    """
    gain = np.ones_like(el_deg)
    taper_mask = el_deg < transition_deg
    if not np.any(taper_mask):
        return gain
    w_trans = np.sin(np.radians(transition_deg))
    w = np.sin(np.radians(el_deg[taper_mask]))
    # w_frac goes from 0 at the horizon to 1 at the transition
    w_frac = np.clip(w / w_trans, 0.0, 1.0)
    # Cosine half-cycle: 0 at w=0 → 1 at w=w_trans
    gain[taper_mask] = 0.5 * (1.0 - np.cos(np.pi * w_frac))
    return gain


@functools.cache
def base_tart_beam(degree=8, q=1, **kw):
    """The base TART antenna beam pattern, pointed at the zenith.

    The beam is **uniform** above ``elevation = 10°`` and **tapers smoothly to
    zero** at the horizon (``elevation = 0°``).  Behind the antenna (below the
    horizon) the response is identically zero.

    The taper is implemented as a cosine half-cycle in ``w = sin(elevation)``,
    giving a smooth C1 transition.  The underlying :class:`Beam` model applies the
    ``w**q`` horizon factor to guarantee a clean zero at the horizon.

    The beam is fitted once and cached by ``(degree, q)`` so repeated calls are
    cheap.

    Parameters
    ----------
    degree : int
        Maximum spherical-harmonic degree for the fit (default 8).
        Higher values capture sharper features but increase the basis size.
    q : int
        Horizon taper power (default 1).  ``q >= 1`` ensures the beam goes
        smoothly (C0 for q=1, C1 for q=2) to zero at the horizon.
    **kw
        Extra keyword arguments forwarded to :class:`Beam` (e.g. ``name``,
        ``frequency_hz``).

    Returns
    -------
    Beam
        A zenith-pointing :class:`Beam` with the base TART primary beam pattern.

    Example
    -------
    >>> from tart_beam.defaults import base_tart_beam
    >>> beam = base_tart_beam(degree=8, q=1, name="TART_dipole")
    >>> import numpy as np
    >>> # Elevations above 10° → uniform response ≈ 1
    >>> zenith = np.array([[0.0, 0.0, 1.0]])
    >>> float(beam.evaluate(zenith))
    1.0
    >>> # On the horizon → zero
    >>> horizon = np.array([[1.0, 0.0, 0.0]])
    >>> float(beam.evaluate(horizon))
    0.0
    """
    s_hat, w, el_deg = _fib_front_hemisphere(n_total=10000)
    gain = _cosine_taper_profile(el_deg, _ELEVATION_TRANSITION_DEG)

    # Name defaults to reflect what this beam represents
    kw.setdefault("name", f"TART_base_N{degree}_q{q}")
    return Beam.fit(s_hat, gain, degree=degree, q=q, **kw)
