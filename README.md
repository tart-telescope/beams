# beams

Beam patterns for TART antennas.

Wide-angle (up to 180° field of view) antenna beams represented as **polynomials
on the sphere** — a truncated spherical-harmonic series in the direction cosines
of a local frame, tapered to zero at the horizon and identically zero behind the
antenna. See [`DESIGN.md`](DESIGN.md) for the full rationale.

## Install

```bash
pip install -e .            # core (numpy, scipy)
pip install -e '.[healpix]' # + healpy for full-sphere combine
```

## Quick start

```python
import numpy as np
from tart_beam import Beam, coverage_map, mosaic, fit_from_json

# fit directly from measured data: a JSON list of {el, az, gain} (degrees),
# with the beam assumed to point at the zenith
beam = fit_from_json("measured.json", degree=8, q=2)

# ...or fit from arbitrary samples (sky unit vectors + values)
beam = Beam.fit(s_hat, values, degree=8, q=2)

# evaluate anywhere; exactly zero behind the antenna
response = beam.evaluate(sky_vectors)

# point copies of it and tile the full HEALPix sphere
beams = [Beam(beam.indices, beam.coeffs, boresight=p, q=beam.q) for p in pointings]
A = coverage_map(beams, nside=64)                 # summed sensitivity
sky, weight = mosaic(beams, per_beam_maps)        # inverse-variance mosaic

beam.to_json("beam.json"); Beam.from_json("beam.json")
```

See [`examples/demo.py`](examples/demo.py) for an end-to-end example.

## Usage

### Default TART beam

The package ships a reference primary-beam pattern for TART antenna elements:
uniform above 10° elevation, tapering smoothly to zero at the horizon.

```python
from tart_beam import base_tart_beam

beam = base_tart_beam(degree=8, q=1)          # first call fits and caches
beam = base_tart_beam()                        # subsequent calls are instant
```

The returned :class:`Beam` is a full spherical-harmonic model, pointable and
compatible with all the operations below.

### Evaluate at elevation / azimuth points

Convert elevation/azimuth (degrees) to unit sky vectors, then call
:meth:`Beam.evaluate`:

```python
from tart_beam import base_tart_beam, elaz_to_vec
import numpy as np

beam = base_tart_beam()

# Single direction
el, az = 45.0, 30.0
s_hat = elaz_to_vec(el, az)           # unit vector on the sphere
gain = beam.evaluate(s_hat[np.newaxis])[0]

# Grid over the whole sky
el = np.linspace(0, 90, 19)          # 0° … 90°
az = np.linspace(0, 360, 37)         # 0° … 360°
EL, AZ = np.meshgrid(el, az, indexing="ij")
s_hat = elaz_to_vec(EL.ravel(), AZ.ravel())
gain = beam.evaluate(s_hat).reshape(EL.shape)
```

### Get a HEALPix map from a beam

Evaluate the beam on a full-sphere HEALPix grid to produce all-sky maps suitable
for mosaicking or visualisation:

```python
from tart_beam import base_tart_beam
from tart_beam import coverage_map, healpix_directions

beam = base_tart_beam()

# Pixel-centre directions for a given resolution (requires healpy)
s_hat = healpix_directions(nside=64)   # shape (49152, 3)
hpx_map = beam.evaluate(s_hat)         # shape (49152,) — the beam on the sphere
```

For overlapping beams pointed in different directions, use the combine helpers
to build coverage maps and mosaics directly on HEALPix:

```python
# Tile four copies of the base beam over the sky
from tart_beam import Beam

pointings = [
    [0.0, 0.0, 1.0],                             # zenith
    [0.866, 0.0, 0.5],                           # 60° from zenith
    [-0.433, 0.75, 0.5],
    [-0.433, -0.75, 0.5],
]
beams = [Beam(beam.indices, beam.coeffs, boresight=p, q=beam.q)
         for p in pointings]

cov = coverage_map(beams, nside=64)              # summed sensitivity
sky, weight = mosaic(beams, per_beam_maps)        # inverse-variance mosaic
```

### Visualise a beam

The ``tart_beam.viz`` module provides helpers for plotting beam patterns
(requires ``matplotlib``; HEALPix plots also need ``healpy``).

```python
from tart_beam import base_tart_beam, viz

beam = base_tart_beam()

# Elevation cut at a fixed azimuth — shows the uniform region and horizon taper
fig, ax = viz.plot_el_cut(beam, az=0.0)

# Azimuth cut at a fixed elevation — reveals rotational symmetry
fig, ax = viz.plot_az_cut(beam, el=45.0)

# Full-sky HEALPix mollweide map (requires healpy)
viz.plot_healpix(beam, nside=64)
```

## Input data format

`fit_from_json` reads a JSON list of records, elevation and azimuth in degrees:

```json
[{"el": 90.0, "az": 0.0, "gain": 1.0},
 {"el": 45.0, "az": 30.0, "gain": 0.42}]
```

The beam is assumed to point at the **zenith**: a sample at elevation `el` maps
to `w = sin(el)`, so `el = 90°` is boresight, `el = 0°` is the horizon, and
samples below the horizon (`el < 0°`) are dropped by the fit. Re-point the
fitted beam afterwards with `beam.set_pointing(boresight)`.

## Layout

- `tart_beam/spherical.py` — real spherical-harmonic basis
- `tart_beam/beam.py` — `Beam`: pointing, evaluation, fitting, (de)serialisation
- `tart_beam/loaders.py` — el/az/gain JSON loader and `fit_from_json`
- `tart_beam/combine.py` — full-sphere HEALPix products (coverage, mosaic, blend)
- `tart_beam/viz.py` — visualisation: elevation/azimuth cuts and HEALPix maps

---
*Copyright (c) 2026 Tim Molteno (tim@elec.ac.nz)*
