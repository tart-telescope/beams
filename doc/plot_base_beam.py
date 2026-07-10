#!/usr/bin/env python3
"""Generate beam visualisations for the documentation.

Saves PNG figures to ``doc/``.  Run from the repo root::

    python doc/plot_base_beam.py
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from tart_beam import base_tart_beam, viz

OUT = "doc"


def _save(fig, name):
    path = f"{OUT}/{name}"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"  saved  {path}")
    plt.close(fig)


def elevation_cut(beam):
    """Beam response vs elevation at azimuth 0°, with annotations."""
    fig, ax = viz.plot_el_cut(beam, az=0.0, linewidth=2.0, color="C0")

    # Annotate the uniform region and taper
    ax.axvspan(10, 90, alpha=0.06, color="C2", label="Uniform region (≥10°)")
    ax.axvspan(0, 10, alpha=0.06, color="C1", label="Horizon taper (0°–10°)")
    ax.axvline(10, color="C2", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.text(12, ax.get_ylim()[1] * 0.92, "10°", color="C2", fontsize=9)

    ax.set_title("Base TART beam — elevation cut at azimuth 0°")
    ax.legend(fontsize=8)
    fig.set_size_inches(6, 4)
    _save(fig, "base_beam_el_cut.png")


def azimuth_cuts(beam):
    """Azimuth cuts at several elevations on the same axes."""
    fig, ax = plt.subplots(figsize=(6, 4))

    for el, ls, label in [
        (90, "-", "Elevation 90° (zenith)"),
        (45, "--", "Elevation 45°"),
        (10, ":", "Elevation 10° (transition)"),
        (5, "-.", "Elevation 5° (in taper)"),
    ]:
        viz.plot_az_cut(beam, el=el, ax=ax, linestyle=ls, linewidth=1.5,
                        label=label)

    ax.set_title("Base TART beam — azimuth cuts")
    ax.legend(fontsize=8)
    fig.tight_layout()
    _save(fig, "base_beam_az_cuts.png")


def healpix_map(beam):
    """Full-sky HEALPix mollview, saved from the healpy figure."""
    try:
        import healpy as hp
    except ImportError:
        print("  (healpy not available — skipping HEALPix map)")
        return

    nside = 64
    from tart_beam.combine import healpix_directions

    s_hat = healpix_directions(nside)
    hpx = beam.evaluate(s_hat)
    hp.mollview(hpx, title="Base TART beam (HEALPix nside=64)", cmap="viridis",
                unit="Response", cbar=True)
    hp.graticule(verbose=False)
    fig = plt.gcf()
    _save(fig, "base_beam_healpix.png")


def main():
    print("Generating base TART beam plots ...")
    # Use a modest degree for clean, smooth plots
    beam = base_tart_beam(degree=8, q=1, name="TART_base")
    print(f"  beam: {beam}")

    elevation_cut(beam)
    azimuth_cuts(beam)
    healpix_map(beam)
    print("Done.")


if __name__ == "__main__":
    main()
