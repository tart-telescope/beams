"""Visualise TART beams on the sky.

Functions for plotting beam patterns on HEALPix maps (via ``healpy``) and
as one-dimensional cuts in elevation or azimuth (via ``matplotlib``).  Both
dependencies are imported lazily so the rest of the package remains lightweight.
"""

# Copyright (c) 2026 Tim Molteno (tim@elec.ac.nz)

import numpy as np

from .loaders import elaz_to_vec


def plot_healpix(beam, nside=64, nest=False, title=None, projection="moll",
                 **hp_kw):
    """Plot the beam on a full-sky HEALPix map.

    Parameters
    ----------
    beam : Beam
        The beam to visualise.
    nside : int
        HEALPix resolution (default 64).
    nest : bool
        Pixel ordering (default False → ring order).
    title : str, optional
        Plot title.
    projection : str
        Healpy projection: ``"moll"`` (mollweide, default), ``"orth"``
        (orthographic), ``"cart"`` (Cartesian), etc.
    **hp_kw
        Extra keyword arguments forwarded to ``healpy.visu.mollview`` (or
        the corresponding projection function).

    Returns
    -------
    Figure
        The healpy figure handle.
    """
    import healpy as hp

    # Evaluate the beam on the HEALPix sphere
    from .combine import healpix_directions

    s_hat = healpix_directions(nside, nest=nest)
    hpx_map = beam.evaluate(s_hat)

    # Build the title from the beam if not given
    if title is None:
        title = repr(beam)

    proj = projection.lower()
    proj_funcs = {
        "moll": hp.mollview,
        "orth": hp.orthview,
        "cart": hp.cartview,
        "gnom": hp.gnomview,
    }
    func = proj_funcs.get(proj, hp.mollview)
    func(hpx_map, title=title, **hp_kw)
    return func


def plot_el_cut(beam, az=0.0, npts=500, ax=None, **plot_kw):
    """Beam response vs elevation at a fixed azimuth.

    Parameters
    ----------
    beam : Beam
    az : float
        Fixed azimuth in degrees (default 0).
    npts : int
        Number of elevation samples from 0° to 90°.
    ax : matplotlib Axes, optional
    **plot_kw
        Forwarded to ``ax.plot``.

    Returns
    -------
    (fig, ax)
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    el = np.linspace(0.0, 90.0, npts)
    s_hat = elaz_to_vec(el, az)
    response = beam.evaluate(s_hat)

    defaults = dict(color="C0", linewidth=1.5)
    defaults.update(plot_kw)
    ax.plot(el, response, **defaults)
    ax.set_xlabel("Elevation (deg)")
    ax.set_ylabel("Beam response")
    ax.set_title(repr(beam))
    ax.set_xlim(0, 90)
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.3)
    return fig, ax


def plot_az_cut(beam, el=45.0, npts=500, ax=None, **plot_kw):
    """Beam response vs azimuth at a fixed elevation.

    Parameters
    ----------
    beam : Beam
    el : float
        Fixed elevation in degrees (default 45).
    npts : int
        Number of azimuth samples from 0° to 360°.
    ax : matplotlib Axes, optional
    **plot_kw
        Forwarded to ``ax.plot``.

    Returns
    -------
    (fig, ax)
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    az = np.linspace(0.0, 360.0, npts)
    s_hat = elaz_to_vec(el, az)
    response = beam.evaluate(s_hat)

    defaults = dict(color="C0", linewidth=1.5)
    defaults.update(plot_kw)
    ax.plot(az, response, **defaults)
    ax.set_xlabel("Azimuth (deg)")
    ax.set_ylabel("Beam response")
    ax.set_title(repr(beam))
    ax.set_xlim(0, 360)
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.3)
    return fig, ax
