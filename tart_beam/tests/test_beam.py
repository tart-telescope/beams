"""Tests for the Beam model: evaluation, horizon behaviour, fitting, pointing."""

import numpy as np
import pytest

from tart_beam import Beam, make_frame, sh_indices


def fib_sphere(n):
    """Roughly uniform unit vectors over the whole sphere."""
    i = np.arange(n) + 0.5
    z = 1.0 - 2.0 * i / n
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, 1.0))
    phi = np.pi * (1.0 + np.sqrt(5.0)) * i
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=-1)


def cos_cubed(s_hat, boresight=(0, 0, 1)):
    """Synthetic axisymmetric beam: w**3 in front, zero behind."""
    w = s_hat @ np.asarray(boresight, float)
    return np.where(w >= 0, np.clip(w, 0, None) ** 3, 0.0)


# -- frame ---------------------------------------------------------------
def test_make_frame_orthonormal():
    for b in ([0, 0, 1], [1, 0, 0], [1, 2, 3], [0.1, 0, 0.99]):
        e1, e2, e3 = make_frame(b)
        E = np.stack([e1, e2, e3])
        assert np.allclose(E @ E.T, np.eye(3), atol=1e-12)
        # e3 points along boresight, right-handed
        assert np.allclose(e3, np.asarray(b, float) / np.linalg.norm(b))
        assert np.allclose(np.cross(e1, e2), e3, atol=1e-12)


# -- evaluation / horizon ------------------------------------------------
def test_zero_behind_antenna():
    beam = Beam([(0, 0)], [1.0], boresight=[0, 0, 1], q=2)
    back = np.array([[0, 0, -1.0], [0.3, 0.2, -0.9], [1, 0, -0.01]])
    back /= np.linalg.norm(back, axis=1, keepdims=True)
    assert np.allclose(beam.evaluate(back), 0.0)


def test_horizon_is_zero_and_continuous():
    beam = Beam([(0, 0)], [3.0], boresight=[0, 0, 1], q=2)
    # exactly on the horizon -> zero because of the w**q taper
    horizon = np.array([[1.0, 0, 0], [0, 1.0, 0], [-1, 0, 0]])
    assert np.allclose(beam.evaluate(horizon), 0.0, atol=1e-12)
    # approaching the horizon from the front -> continuous decay to zero
    eps = np.array([0.1, 0.01, 0.001])
    s = np.stack([np.sqrt(1 - eps ** 2), np.zeros_like(eps), eps], axis=-1)
    vals = beam.evaluate(s)
    assert np.all(vals > 0)
    assert vals[0] > vals[1] > vals[2]


def test_q_controls_taper():
    a, b = (Beam([(0, 0)], [1.0], q=q) for q in (1, 3))
    s = np.array([[np.sqrt(1 - 0.04), 0, 0.2]])  # w = 0.2
    # higher q suppresses the response harder near the horizon
    assert b.evaluate(s)[0] < a.evaluate(s)[0]


def test_evaluate_shape_preserved():
    beam = Beam([(0, 0)], [1.0])
    s = fib_sphere(50).reshape(10, 5, 3)
    assert beam.evaluate(s).shape == (10, 5)


# -- fitting -------------------------------------------------------------
def test_fit_recovers_axisymmetric_beam():
    s = fib_sphere(8000)
    vals = cos_cubed(s)
    beam = Beam.fit(s, vals, degree=8, q=2, ridge=1e-8)
    assert beam.residual_rms(s, vals) < 1e-3
    # only m == 0 harmonics should carry weight (axisymmetric truth)
    big = [(l, m) for (l, m), a in zip(beam.indices, beam.coeffs)
           if abs(a) > 1e-2]
    assert all(m == 0 for (l, m) in big)


def test_fit_then_zero_behind():
    s = fib_sphere(4000)
    beam = Beam.fit(s, cos_cubed(s), degree=6, q=2)
    back = np.array([[0, 0, -1.0], [0.5, 0.5, -0.707]])
    back /= np.linalg.norm(back, axis=1, keepdims=True)
    assert np.allclose(beam.evaluate(back), 0.0)


def test_weights_accepted():
    s = fib_sphere(2000)
    vals = cos_cubed(s)
    w = np.ones(len(s))
    beam = Beam.fit(s, vals, degree=5, q=2, weights=w)
    assert beam.residual_rms(s, vals) < 1e-2


# -- pointing invariance -------------------------------------------------
def test_pointing_invariance():
    """A pointed beam equals the canonical beam evaluated in rotated coords."""
    s = fib_sphere(3000)
    canon = Beam.fit(s, cos_cubed(s), degree=6, q=2)

    p = np.array([0.3, -0.4, 0.866])
    p /= np.linalg.norm(p)
    pointed = Beam(canon.indices, canon.coeffs, boresight=p, q=canon.q)

    # response toward the new boresight should match response toward old one
    assert np.isclose(pointed.evaluate(p[None])[0],
                      canon.evaluate(np.array([[0, 0, 1.0]]))[0], atol=1e-6)
    # and the pointed beam is zero opposite its boresight
    assert np.isclose(pointed.evaluate((-p)[None])[0], 0.0)


# -- serialisation -------------------------------------------------------
def test_json_round_trip(tmp_path):
    s = fib_sphere(2000)
    beam = Beam.fit(s, cos_cubed(s), degree=5, q=2,
                    name="rt", frequency_hz=1.5e9)
    path = tmp_path / "beam.json"
    beam.to_json(path)
    other = Beam.from_json(path)

    assert other.name == "rt"
    assert other.frequency_hz == 1.5e9
    assert other.q == beam.q
    assert other.indices == beam.indices
    assert np.allclose(other.coeffs, beam.coeffs)
    assert np.allclose(other.evaluate(s), beam.evaluate(s))


def test_length_mismatch_raises():
    with pytest.raises(ValueError):
        Beam([(0, 0), (1, 0)], [1.0])
