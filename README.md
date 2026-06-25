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
