"""tart_beam: wide-angle polynomial beam patterns on the sphere.

See ``DESIGN.md`` for the representation: a truncated spherical-harmonic series
in the direction cosines of a local frame, tapered by ``w**q`` and zero behind
the antenna.
"""

from .beam import Beam, make_frame
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

__version__ = "0.1.0"

__all__ = [
    "Beam",
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
]
