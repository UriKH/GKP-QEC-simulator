"""Optimal (entanglement-fidelity-maximizing) recovery for a qubit code: an SDP, solved with cvxpy/Clarabel.

Why: Petz is near-optimal (Barnum-Knill: within a factor 2 in entanglement infidelity) but not optimal; for the
bare Fock qubit it is even beaten by plain projection. To compare codes at their best, we also compute the
optimal recovery R for the given channel E and encoder V.

F_e(R) = Tr[J_R X] is linear in the Choi matrix J_R = sum_ab R(|a><b|) (x) |a><b|  (2N x 2N, output (x) input),
with X = (1/4) sum_ij |j><i| (x) sigma_ij^T,  sigma_ij = E(V|i><j|V^dag).  The constraint is Tr_out J_R = 1_N.
Maximizing a linear function over Choi matrices is a semidefinite program (Fletcher, Shor & Win, PRA 75, 012338
(2007); Kosut & Lidar 2009). We solve it with cvxpy + Clarabel restricted to the input support S of X (the span of
everything the channel can output; the neglected weight Tr[X (1 - P_S)] bounds the error), and certify the result
with an independent dual-feasible point (`upper_bound`, weak duality).
`optimize` (the Reimpell-Werner iteration, PRL 94, 080501 (2005)) is kept for reference only: in our tests it
stalled at start-dependent values, so it is NOT used for results.

cvxpy is not in the qcsim environment. It is installed outside it (pip --target <dir>) and found via the
CVXPY_PATH environment variable; see results_ext/README.md.

This lives outside cvqec/ on purpose (extension work); it reuses the library's channels and codes unchanged.
"""
import numpy as np


def choi_X(code, channel, N):
    V = code.encoder(N).full()
    import qutip as qt
    X = np.zeros((2 * N, 2 * N), complex)
    sig = {}
    for i in range(2):
        for j in range(2):
            s = channel.apply(qt.Qobj(np.outer(V[:, i], V[:, j].conj()))).full()
            sig[i, j] = s
            Eji = np.zeros((2, 2)); Eji[j, i] = 1
            X += np.kron(Eji, s.T) / 4
    return (X + X.conj().T) / 2, sig


def _ptrace_out(J, N):
    return J.reshape(2, N, 2, N).trace(axis1=0, axis2=2)


def _inv_sqrt(L, rcond=1e-14):
    w, U = np.linalg.eigh((L + L.conj().T) / 2)
    keep = w > rcond * w.max()
    return (U[:, keep] / np.sqrt(w[keep])) @ U[:, keep].conj().T


def optimize(X, N, iters=2000, tol=1e-12, J0=None):
    J = np.kron(np.eye(2) / 2, np.eye(N)) if J0 is None else J0
    f_old = np.real(np.trace(J @ X))
    hist = [f_old]
    for k in range(iters):
        Y = J @ X @ J
        S = np.kron(np.eye(2), _inv_sqrt(_ptrace_out(Y, N)))
        J = S @ Y @ S
        J = (J + J.conj().T) / 2
        f = np.real(np.trace(J @ X))
        hist.append(f)
        if abs(f - f_old) < tol:
            break
        f_old = f
    return J, f, hist


def upper_bound(J, X, N):
    """Certified upper bound on max F_e from the dual SDP:  min Tr[Y]  s.t.  1_2 (x) Y >= X.

    Candidate Y = Hermitian part of Tr_out(X J) (exact at the optimum, by the KKT conditions). If
    M = 1 (x) Y - X has negative eigenvalues -l_k (eigenvectors v_k), lift Y by 2 sum_k l_k Tr_out|v_k><v_k|:
    since rho_AB <= d_A (1 (x) rho_B) for any state (the map X -> d Tr(X) 1 - X is completely positive),
    the lifted Y is dual feasible, and Tr[Y] grows only by 2 sum_k l_k. Weak duality: F_e(opt) <= Tr[Y].
    """
    Y = _ptrace_out(X @ J, N)
    Y = (Y + Y.conj().T) / 2
    w, U = np.linalg.eigh(np.kron(np.eye(2), Y) - X)
    neg = w < 0
    lift = sum(-wk * _ptrace_out(np.outer(U[:, k], U[:, k].conj()), N) for k, wk in zip(np.where(neg)[0], w[neg]))
    Yf = Y + 2 * lift if neg.any() else Y
    return float(np.real(np.trace(Yf))), float(-w[neg].sum()) if neg.any() else 0.0


def best_recovery(code, channel, N):
    """Optimal recovery for (code, channel): SDP optimum and a certified upper bound."""
    X, _ = choi_X(code, channel, N)
    return solve_sdp(X, N)


def solve_sdp(X, N, support_tol=1e-13):
    """Optimal F_e by SDP on out (x) S. Returns dict with the optimum, the certified dual bound and diagnostics."""
    import os, sys
    if os.environ.get("CVXPY_PATH") and os.environ["CVXPY_PATH"] not in sys.path:
        sys.path.append(os.environ["CVXPY_PATH"])     # appended: numpy/scipy still come from qcsim
    import cvxpy as cp
    w, U = np.linalg.eigh(_ptrace_out(X, N))
    keep = w > support_tol * w.max()
    B = U[:, keep]                                   # N x d basis of S
    d = B.shape[1]
    T = np.kron(np.eye(2), B)                        # (2N) x (2d)
    Xs = T.conj().T @ X @ T
    Xs = (Xs + Xs.conj().T) / 2
    neglected = float(np.real(np.trace(X)) - np.real(np.trace(Xs)))
    J = cp.Variable((2 * d, 2 * d), hermitian=True)
    cons = [J >> 0, cp.partial_trace(J, (2, d), axis=0) == np.eye(d)]
    prob = cp.Problem(cp.Maximize(cp.real(cp.trace(J @ Xs))), cons)
    prob.solve(solver=cp.CLARABEL)
    Jv = J.value
    Jv = (Jv + Jv.conj().T) / 2
    # repair tiny infeasibility (solver tolerance): clip negative eigenvalues, re-normalize the input marginal
    ev, EV = np.linalg.eigh(Jv)
    Jv = (EV * np.clip(ev, 0, None)) @ EV.conj().T
    S = np.kron(np.eye(2), _inv_sqrt(_ptrace_out(Jv, d)))
    Jv = S @ Jv @ S
    f = float(np.real(np.trace(Jv @ Xs)))           # achieved by an exactly valid recovery (lower bound)
    ub, _ = upper_bound(Jv, Xs, d)
    # recovery on all of N: J_full = T Jv T^dag + (1/2) (x) (1 - P_S); the neglected part can add at most `neglected`
    return {"F_e_opt": f, "F_e_opt_upper": ub + neglected, "support_dim": int(d), "neglected_weight": neglected,
            "solver_status": prob.status, "F_avg_opt": (2 * f + 1) / 3}
