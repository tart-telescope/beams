"""tart-beam: wide-angle polynomial beam patterns on the sphere.

See ``DESIGN.md`` for the representation: a truncated spherical-harmonic series
in the direction cosines of a local frame, tapered by ``w**q`` and zero behind
the antenna.
"""

# Copyright (c) 2026 Tim Molteno (tim@elec.ac.nz)

from .beam import Beam, make_frame
from .defaults import base_tart_beam
from .combine import (
    coverage_map,
    healpix_directions,
    mosaic,
    partition_of_unity,
)
from .loaders import (
    elaz_to_vec,
    fit_from_json,
    load_elaz_json,
    samples_from_records,
)
from .spherical import basis_matrix, real_sph_harm, sh_indices
from . import viz  # optional (requires healpy + matplotlib)

__version__ = "0.2.0"

__all__ = [
    "Beam",
    "base_tart_beam",
    "make_frame",
    "coverage_map",
    "partition_of_unity",
    "mosaic",
    "healpix_directions",
    "real_sph_harm",
    "basis_matrix",
    "sh_indices",
    "elaz_to_vec",
    "samples_from_records",
    "load_elaz_json",
    "fit_from_json",
    "viz",
]
