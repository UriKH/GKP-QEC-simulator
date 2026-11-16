"""Error-subspace recovery: the Knill-Laflamme picture (course p.74-77) turned into a decoder.

For a list of correctable errors E_k, the damaged codewords E_k|mu_L> span "error subspaces".
Recovery measures which subspace the state is in (for the cat this is the photon-number parity
measurement) and rotates that subspace back onto the code:   W_k = sum_mu |mu><e_{k,mu}|.
Anything outside every error subspace is an uncorrectable error; it is mapped to the maximally
mixed qubit, so the decoder is trace preserving and failures count fully against fidelity.
"""

import numpy as np
import qutip as qt

from cvqec.decoders.base import Decoder, register_decoder


def _loewdin(A: np.ndarray) -> np.ndarray:
    eigval, eigvec = np.linalg.eigh(A.conj().T @ A)
    return A @ eigvec @ np.diag(eigval ** -0.5) @ eigvec.conj().T


class ErrorSubspaceDecoder(Decoder):
    name = "error_subspace"

    def __init__(self, code, errors=None):
        super().__init__(code)
        self.errors = errors or [lambda N: qt.qeye(N)]
        self._cache = {}

    def _isometries(self, N):
        if N not in self._cache:
            V = self.code.encoder(N).full()
            taken = np.zeros((N, 0))
            Ws = []
            for E in self.errors:
                B = E(N).full() @ V
                B = B - taken @ (taken.conj().T @ B)          # keep error subspaces orthogonal
                if np.linalg.norm(B) < 1e-10:
                    continue
                B = _loewdin(B)
                Ws.append(B.conj().T)
                taken = np.hstack([taken, B])
            Q = np.eye(N) - taken @ taken.conj().T
            self._cache[N] = (Ws, Q)
        return self._cache[N]

    def decode(self, rho):
        r = rho.full()
        Ws, Q = self._isometries(r.shape[0])
        out = sum(W @ r @ W.conj().T for W in Ws)
        return out + np.trace(Q @ r) * np.eye(2) / 2


@register_decoder("projection")
class ProjectionDecoder(ErrorSubspaceDecoder):
    """No correction at all: keep the code-space part. The baseline for 'does correction help?'."""

    name = "projection"

    def __init__(self, code):
        super().__init__(code)


@register_decoder("cat_parity")
class CatParityDecoder(ErrorSubspaceDecoder):
    """4-leg cat: errors {1, a}. The two error subspaces have opposite photon-number parity (p.77)."""

    name = "cat_parity"

    def __init__(self, code):
        super().__init__(code, errors=[lambda N: qt.qeye(N), lambda N: qt.destroy(N)])
