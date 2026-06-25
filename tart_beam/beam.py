"""Pointable, polynomial wide-angle beam on the sphere.

A :class:`Beam` is a truncated spherical-harmonic series in the direction
cosines of a local frame, multiplied by a ``w**q`` horizon taper and clipped to
zero behind the antenna (``w < 0``). See ``DESIGN.md``.
"""

import json

import numpy as np

from .spherical import basis_matrix, sh_indices


def make_frame(boresight, roll_ref=None):
    """Build an orthonormal frame ``(e1, e2, e3)`` with ``e3`` along boresight.

    ``roll_ref`` fixes the azimuth origin (``e1`` direction). It defaults to the
    global x-axis, or y-axis when boresight is too close to x.
    """
    e3 = np.asarray(boresight, dtype=float)
    e3 = e3 / np.linalg.norm(e3)
    if roll_ref is None:
        roll_ref = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(roll_ref, e3)) > 0.9:
            roll_ref = np.array([0.0, 1.0, 0.0])
    roll_ref = np.asarray(roll_ref, dtype=float)
    e1 = roll_ref - np.dot(roll_ref, e3) * e3
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(e3, e1)
    return e1, e2, e3


class Beam:
    """A wide-angle beam pattern pointed at ``boresight`` on the sphere.

    Parameters
    ----------
    indices : list[(l, m)]
        Spherical-harmonic indices, parallel to ``coeffs``.
    coeffs : array
        Real spherical-harmonic coefficients ``a_l^m``.
    boresight, roll_ref : array-like
        Define the local frame (see :func:`make_frame`).
    q : int
        Horizon taper power (``B = w**q * core``); ``q >= 1`` forces a smooth
        zero at the horizon, ``q >= 2`` makes it C1.
    """

    def __init__(self, indices, coeffs, boresight=(0.0, 0.0, 1.0), roll_ref=None,
                 q=2, name=None, frequency_hz=None):
        self.indices = [tuple(i) for i in indices]
        self.coeffs = np.asarray(coeffs, dtype=float)
        if len(self.indices) != self.coeffs.size:
            raise ValueError("indices and coeffs must have equal length")
        self.q = int(q)
        self.name = name
        self.frequency_hz = frequency_hz
        self.set_pointing(boresight, roll_ref)

    # -- pointing ---------------------------------------------------------
    def set_pointing(self, boresight, roll_ref=None):
        """Re-point the beam. Coefficients are unchanged (pointing-invariant)."""
        self.e1, self.e2, self.e3 = make_frame(boresight, roll_ref)
        return self

    @property
    def boresight(self):
        return self.e3

    def local_cosines(self, s_hat):
        """Direction cosines ``(u, v, w)`` of sky vectors in the beam frame."""
        s_hat = np.asarray(s_hat, dtype=float)
        u = s_hat @ self.e1
        v = s_hat @ self.e2
        w = s_hat @ self.e3
        return u, v, w

    # -- evaluation -------------------------------------------------------
    def evaluate(self, s_hat):
        """Beam response at unit sky vectors ``s_hat`` (shape ``(..., 3)``).

        Returns zero exactly where ``w < 0`` (behind the antenna).
        """
        s_hat = np.asarray(s_hat, dtype=float)
        u, v, w = self.local_cosines(s_hat)
        front = w >= 0.0
        # angles measured from boresight
        theta = np.arccos(np.clip(w, -1.0, 1.0))
        phi = np.arctan2(v, u)
        core = basis_matrix(self.indices, theta, phi) @ self.coeffs
        core = core.reshape(w.shape)
        taper = np.where(front, np.clip(w, 0.0, None) ** self.q, 0.0)
        return np.where(front, taper * core, 0.0)

    # -- fitting ----------------------------------------------------------
    @classmethod
    def fit(cls, s_hat, values, degree, q=2, ridge=0.0, weights=None,
            boresight=(0.0, 0.0, 1.0), roll_ref=None, **kw):
        """Least-squares fit a beam to samples ``values`` at directions ``s_hat``.

        Samples behind the antenna (``w < 0``) are dropped. The ``w**q`` taper is
        folded into the design matrix, so no division is needed near the horizon.

        Parameters
        ----------
        degree : int
            Maximum spherical-harmonic degree ``N``.
        ridge : float
            Tikhonov regularisation strength ``lambda``.
        weights : array, optional
            Per-sample weights (e.g. inverse variance).
        """
        indices = sh_indices(degree)
        proto = cls(indices, np.zeros(len(indices)), boresight, roll_ref, q=q)
        u, v, w = proto.local_cosines(s_hat)
        values = np.asarray(values, dtype=float).ravel()

        front = w >= 0.0
        theta = np.arccos(np.clip(w[front], -1.0, 1.0))
        phi = np.arctan2(v[front], u[front])
        # model is  B = w**q * (Y @ a)  ->  fold taper into the design matrix
        M = basis_matrix(indices, theta, phi) * (w[front] ** q)[:, None]
        b = values[front]

        if weights is not None:
            sw = np.sqrt(np.asarray(weights, dtype=float).ravel()[front])
            M = M * sw[:, None]
            b = b * sw

        if ridge > 0.0:
            A = M.T @ M + ridge * np.eye(M.shape[1])
            coeffs = np.linalg.solve(A, M.T @ b)
        else:
            coeffs, *_ = np.linalg.lstsq(M, b, rcond=None)

        return cls(indices, coeffs, boresight, roll_ref, q=q, **kw)

    def residual_rms(self, s_hat, values):
        """RMS of (model - values) over front-hemisphere samples."""
        _, _, w = self.local_cosines(s_hat)
        front = w >= 0.0
        pred = self.evaluate(s_hat)
        resid = pred.ravel()[front] - np.asarray(values, float).ravel()[front]
        return float(np.sqrt(np.mean(resid ** 2)))

    # -- serialisation ----------------------------------------------------
    def to_dict(self):
        return {
            "name": self.name,
            "frequency_hz": self.frequency_hz,
            "basis": "real_spherical_harmonic",
            "degree_N": max(l for l, _ in self.indices) if self.indices else 0,
            "horizon_power_q": self.q,
            "boresight": self.e3.tolist(),
            "roll_axis_e1": self.e1.tolist(),
            "coeffs": [
                {"l": l, "m": m, "a": float(a)}
                for (l, m), a in zip(self.indices, self.coeffs)
            ],
        }

    def to_json(self, path):
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def from_dict(cls, d):
        indices = [(c["l"], c["m"]) for c in d["coeffs"]]
        coeffs = [c["a"] for c in d["coeffs"]]
        return cls(indices, coeffs, boresight=d["boresight"],
                   roll_ref=d.get("roll_axis_e1"), q=d.get("horizon_power_q", 2),
                   name=d.get("name"), frequency_hz=d.get("frequency_hz"))

    @classmethod
    def from_json(cls, path):
        with open(path) as f:
            return cls.from_dict(json.load(f))

    def __repr__(self):
        N = max(l for l, _ in self.indices) if self.indices else 0
        return (f"Beam(name={self.name!r}, N={N}, q={self.q}, "
                f"boresight={np.round(self.e3, 3).tolist()})")
