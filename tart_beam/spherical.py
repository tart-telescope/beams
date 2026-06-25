"""Real spherical-harmonic basis on the sphere.

A beam is a polynomial in the direction cosines ``(u, v, w)`` of a local frame.
Restricted to the unit sphere that polynomial is a truncated spherical-harmonic
series, so we use real spherical harmonics ``Y_l^m(theta, phi)`` as the basis.
See ``DESIGN.md`` section 3.
"""

import numpy as np
from scipy.special import lpmv


def sh_indices(degree):
    """Ordered list of ``(l, m)`` pairs for all harmonics up to ``degree``."""
    return [(l, m) for l in range(degree + 1) for m in range(-l, l + 1)]


def _norm(l, m):
    """Normalisation constant for the real spherical harmonic ``Y_l^m``."""
    from scipy.special import gammaln

    am = abs(m)
    # sqrt((2l+1)/(4 pi) * (l-|m|)! / (l+|m|)!), computed via log-gamma for safety
    log_ratio = gammaln(l - am + 1) - gammaln(l + am + 1)
    return np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(log_ratio))


def real_sph_harm(l, m, theta, phi):
    """Evaluate the real spherical harmonic ``Y_l^m`` at ``(theta, phi)``.

    ``theta`` is the polar angle from boresight, ``phi`` the azimuth. Uses the
    Condon-Shortley convention consistently (via :func:`scipy.special.lpmv`), so
    the basis is internally consistent for fitting and evaluation.
    """
    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)
    am = abs(m)
    leg = lpmv(am, l, np.cos(theta))
    norm = _norm(l, m)
    if m > 0:
        return np.sqrt(2.0) * norm * leg * np.cos(am * phi)
    if m < 0:
        return np.sqrt(2.0) * norm * leg * np.sin(am * phi)
    return norm * leg


def basis_matrix(indices, theta, phi):
    """Design matrix ``M[i, j] = Y_{l_j}^{m_j}(theta_i, phi_i)``.

    Shape ``(n_points, n_basis)``.
    """
    theta = np.asarray(theta, dtype=float).ravel()
    phi = np.asarray(phi, dtype=float).ravel()
    cols = [real_sph_harm(l, m, theta, phi) for (l, m) in indices]
    return np.stack(cols, axis=-1)
