# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

- **`tart_beam/defaults.py`** — new module providing `base_tart_beam()`, the
  nominal TART primary beam pattern. The beam is uniform above 10° elevation
  and tapers smoothly to zero at the horizon via a cosine half-cycle in
  `w = sin(elevation)`. The beam is fitted once and cached by `(degree, q)`
  so repeated calls return the pre-computed result instantly.
  (Issue #… / PR #…)

- **`README.md` documentation** — new "Usage" section covering:
  - Creating a beam with `base_tart_beam()`
  - Evaluating the beam at elevation / azimuth points via `elaz_to_vec`
  - Producing HEALPix all-sky maps from a beam via `healpix_directions`
  - Tiling multiple pointed beams and combining with `coverage_map` / `mosaic`

- **`.github/workflows/publish.yml`** — CI workflow for publishing to PyPI
  via trusted (OIDC) publishing. Triggered by tags matching `v*.*.*`. Runs
  tests, builds sdist + wheel with hatchling, and publishes with
  `pypa/gh-action-pypi-publish@release/v1`.

- **Copyright headers** — all source files now carry
  `Copyright (c) 2026 Tim Molteno (tim@elec.ac.nz)`.

### Changed

- **Package name** — `pyproject.toml` project name changed from `tart_beam`
  to `tart-beam` (the PyPI distribution name). The Python import path
  remains `tart_beam`.
