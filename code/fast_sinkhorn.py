"""Batched Sinkhorn replacement for POT's ``ot.sinkhorn2``.

POT spends about 20 ms of pure Python overhead on every call, which makes a
single Ricci-flow step on a 37-node graph cost several seconds even though the
transport problems themselves are 6x6. This module performs exactly the same
Sinkhorn-Knopp iteration as ``ot.bregman.sinkhorn_knopp`` (and hence the same
loss as ``ot.sinkhorn2``), but batches all edges that share the same support
shape into one numpy einsum pipeline.

Reference algorithm (POT 0.9.7, ot.bregman.sinkhorn_knopp)::

    K  = exp(-M / reg)
    Kp = (1/a)[:, None] * K
    u  = ones(dim_a) / dim_a ; v = ones(dim_b) / dim_b
    repeat:
        KtransposeU = K.T @ u
        v = b / KtransposeU
        u = 1 / (Kp @ v)
        if any(KtransposeU == 0) or any non-finite u/v:
            revert to the previous iterate and stop
        every 10th iteration:
            err = ||einsum('i,ij,j->j', u, K, v) - b||
            if err < stopThr: stop
    G    = u[:, None] * K * v[None, :]
    loss = sum(M * G)

The batched version applies the same update to every edge of a group and freezes
an edge as soon as its own convergence test passes, so the result is identical
to running the reference loop edge by edge.
"""

from __future__ import print_function

import numpy as np

MAX_GROUP_CELLS = 20000000


def _solve_group(a, b, c, reg, num_iter_max, stop_thr):
    batch, dim_a, dim_b = c.shape
    k = np.exp(-c / float(reg))
    kp = k / a[:, :, None]

    u = np.full((batch, dim_a), 1.0 / dim_a, dtype=float)
    v = np.full((batch, dim_b), 1.0 / dim_b, dtype=float)
    active = np.ones(batch, dtype=bool)

    for ii in range(int(num_iter_max)):
        if not active.any():
            break
        k_transpose_u = np.einsum("bij,bi->bj", k, u)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            v_new = b / k_transpose_u
            u_new = 1.0 / np.einsum("bij,bj->bi", kp, v_new)

        zero_hit = np.any(k_transpose_u == 0.0, axis=1)
        finite = (np.all(np.isfinite(v_new), axis=1)
                  & np.all(np.isfinite(u_new), axis=1))
        good = finite & ~zero_hit

        update = active & good
        u = np.where(update[:, None], u_new, u)
        v = np.where(update[:, None], v_new, v)
        active = active & good

        if ii % 10 == 0:
            marginal = np.einsum("bi,bij,bj->bj", u, k, v)
            err = np.linalg.norm(marginal - b, axis=1)
            active = active & ~(err < stop_thr)

    return np.einsum("bi,bij,bj,bij->b", u, k, v, c)


def batched_sinkhorn_cost(a_list, b_list, c_list, reg, num_iter_max=5000,
                          stop_thr=1e-9):
    """Return the entropic transport loss for every edge transport problem."""
    count = len(a_list)
    out = np.zeros(count, dtype=float)
    groups = {}
    for index in range(count):
        shape = (int(a_list[index].shape[0]), int(b_list[index].shape[0]))
        groups.setdefault(shape, []).append(index)

    for shape, indices in groups.items():
        dim_a, dim_b = shape
        chunk = max(1, int(MAX_GROUP_CELLS // max(1, dim_a * dim_b)))
        for start in range(0, len(indices), chunk):
            part = indices[start:start + chunk]
            a = np.asarray([a_list[i] for i in part], dtype=float)
            b = np.asarray([b_list[i] for i in part], dtype=float)
            c = np.asarray([c_list[i] for i in part], dtype=float)
            out[part] = _solve_group(a, b, c, reg, num_iter_max, stop_thr)
    return out
