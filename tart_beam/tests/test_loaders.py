"""Tests for the el/az/gain JSON loader."""

import json

import numpy as np
import pytest

from tart_beam import elaz_to_vec, fit_from_json, load_elaz_json, samples_from_records


def test_elaz_to_vec_cardinal_points():
    # zenith
    assert np.allclose(elaz_to_vec(90.0, 0.0), [0, 0, 1], atol=1e-12)
    # horizon along +x (az=0) and +y (az=90)
    assert np.allclose(elaz_to_vec(0.0, 0.0), [1, 0, 0], atol=1e-12)
    assert np.allclose(elaz_to_vec(0.0, 90.0), [0, 1, 0], atol=1e-12)
    # below horizon -> negative z (behind a zenith-pointing antenna)
    assert elaz_to_vec(-30.0, 0.0)[2] < 0


def test_elaz_to_vec_unit_norm():
    el = np.array([90, 45, 10, 0, -20.0])
    az = np.array([0, 30, 200, 359, 12.0])
    v = elaz_to_vec(el, az)
    assert np.allclose(np.linalg.norm(v, axis=-1), 1.0)


def test_samples_from_records():
    recs = [{"el": 90, "az": 0, "gain": 1.0}, {"el": 0, "az": 0, "gain": 0.0}]
    s, g = samples_from_records(recs)
    assert s.shape == (2, 3)
    assert np.allclose(g, [1.0, 0.0])


def test_samples_from_records_empty_raises():
    with pytest.raises(ValueError):
        samples_from_records([])


def _synthetic_records(n=2000):
    """A cos^3 zenith beam sampled on an el/az grid, including below horizon."""
    rng = np.random.default_rng(0)
    el = rng.uniform(-90, 90, n)
    az = rng.uniform(0, 360, n)
    w = np.sin(np.radians(el))            # = cos(theta from zenith)
    gain = np.where(w >= 0, np.clip(w, 0, None) ** 3, 0.0)
    return [{"el": float(e), "az": float(a), "gain": float(g)}
            for e, a, g in zip(el, az, gain)]


def test_fit_from_json_round_trip(tmp_path):
    path = tmp_path / "beam_data.json"
    path.write_text(json.dumps(_synthetic_records()))

    s_hat, gain = load_elaz_json(path)
    assert s_hat.shape[0] == gain.shape[0]

    beam = fit_from_json(path, degree=8, q=2, name="loaded")
    # boresight is the zenith, as documented
    assert np.allclose(beam.boresight, [0, 0, 1])
    assert beam.name == "loaded"
    assert beam.residual_rms(s_hat, gain) < 1e-2
    # exactly zero behind the antenna
    assert np.isclose(beam.evaluate(np.array([[0, 0, -1.0]]))[0], 0.0)
