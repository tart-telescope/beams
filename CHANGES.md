# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

### Changed

### Fixed

---

## [0.2.0] — 2026-07-10

### Added

- **`tart_beam/viz.py`** — new module for beam visualisation: `plot_healpix()`
  (full-sky HEALPix map via healpy), `plot_el_cut()` (beam response vs elevation
  at a fixed azimuth), and `plot_az_cut()` (beam response vs azimuth at a fixed
  elevation).  Both ``matplotlib`` and ``healpy`` are imported lazily so the core
  package remains lightweight.

- **`doc/` directory** — script `doc/plot_base_beam.py` and generated PNG figures
  (elevation cut, azimuth cuts, HEALPix map).  Referenced from the README in a
  visualisation gallery table.

- **README.md visualisation docs** — new subsection documenting `tart_beam.viz`
  usage with embedded example plots.

### Changed

- **`testing.yaml` removed** — unit tests are now run only as part of the
  `publish.yml` workflow, triggered by `v*.*.*` tags.

### Fixed

- **`elaz_to_vec` broadcasting** — fixed a crash when passing a scalar elevation
  with an array of azimuths (or vice versa).  Now uses `np.broadcast_arrays` to
  handle mixed scalar/array inputs.

---

## [0.1.0] — 2026-07-10

### Added

- **`tart_beam/defaults.py`** — new module providing `base_tart_beam()`, the
  nominal TART primary beam pattern. The beam is uniform above 10° elevation
  and tapers smoothly to zero at the horizon via a cosine half-cycle in
  `w = sin(elevation)`. The beam is fitted once and cached by `(degree, q)`
  so repeated calls return the pre-computed result instantly.

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
