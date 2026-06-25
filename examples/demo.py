"""Demo: fit a wide-angle beam, then tile the sphere with overlapping beams.

Run from the repo root:  python examples/demo.py
The HEALPix steps are skipped automatically if healpy is not installed.
"""

import os

import numpy as np

from tart_beam import Beam, fit_from_json, load_elaz_json, sh_indices

HERE = os.path.dirname(__file__)
SAMPLE_JSON = os.path.join(HERE, "sample_beam.json")


def sample_directions(n):
    """Roughly uniform unit vectors over the whole sphere (Fibonacci sphere)."""
    i = np.arange(n) + 0.5
    z = 1.0 - 2.0 * i / n
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, 1.0))
    phi = np.pi * (1.0 + np.sqrt(5.0)) * i
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=-1)


def truth(s_hat, boresight=(0, 0, 1)):
    """A synthetic 'measured' beam: cos^3 fall-off with a mild azimuth lobe."""
    w = s_hat @ np.asarray(boresight, float)
    u = s_hat @ np.array([1.0, 0.0, 0.0])
    front = w >= 0.0
    val = (np.clip(w, 0, None) ** 3) * (1.0 + 0.3 * (u ** 2))
    return np.where(front, val, 0.0)


def main():
    # --- 0. fit straight from an el/az/gain JSON file (beam at zenith) -------
    print(f"Loading measured samples from {os.path.basename(SAMPLE_JSON)}")
    s_meas, g_meas = load_elaz_json(SAMPLE_JSON)
    loaded = fit_from_json(SAMPLE_JSON, degree=6, q=2, name="from_json")
    print(f"  {loaded}  (residual RMS {loaded.residual_rms(s_meas, g_meas):.3e})\n")

    # --- 1. fit a beam to synthetic samples ---------------------------------
    s = sample_directions(6000)
    measured = truth(s)

    beam = Beam.fit(s, measured, degree=8, q=2, ridge=1e-6,
                    name="demo", frequency_hz=1.57542e9)
    print(beam)
    print(f"  basis size: {len(sh_indices(8))} harmonics")
    print(f"  fit residual RMS: {beam.residual_rms(s, measured):.3e}")

    # exactly zero behind the antenna
    back = np.array([[0.0, 0.0, -1.0], [0.2, 0.1, -0.9]])
    print(f"  response behind antenna: {beam.evaluate(back)}  (must be 0)")

    # round-trip through JSON
    beam.to_json("/tmp/demo_beam.json")
    reloaded = Beam.from_json("/tmp/demo_beam.json")
    print(f"  reloaded: {reloaded}")

    # --- 2. point copies of the beam in different directions ----------------
    def unit(theta_deg, phi_deg):
        t, p = np.radians(theta_deg), np.radians(phi_deg)
        return [np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)]

    pointings = [unit(0, 0), unit(60, 0), unit(60, 120), unit(60, 240)]
    beams = [Beam(beam.indices, beam.coeffs, boresight=p, q=beam.q)
             for p in pointings]

    # --- 3. combine on the full HEALPix sphere ------------------------------
    try:
        from tart_beam import coverage_map, partition_of_unity
        nside = 32
        cov = coverage_map(beams, nside)
        weights = partition_of_unity(beams, nside)
        covered = (cov > 0).mean()
        print(f"\nHEALPix nside={nside}: sky fraction covered = {covered:.2%}")
        print(f"  partition-of-unity weights sum (covered px) ~ "
              f"{weights.sum(axis=0)[cov > 0].mean():.4f}")
    except ImportError:
        print("\n(healpy not installed - skipping full-sphere combine)")


if __name__ == "__main__":
    main()
