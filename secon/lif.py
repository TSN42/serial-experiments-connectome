"""Shiu et al. (2024) 型の LIF ネットワークを numba で回す。

    dv/dt = (v0 - v + g - a) / tau_m  （不応期中は v, g とも固定 = Brian2 の "unless refractory"）
    dg/dt = -g / tau_syn
    da/dt = -a / tau_adapt、発火ごとに a += adapt_b（スパイク頻度順応。adapt_b = 0 なら Shiu 2024 と同じ）
    v > v_th で発火 → v = v_rst、t_ref の不応期。シナプス後の g へは t_delay 後に w_syn × 符号付きシナプス数を加算。
    刺激ニューロンには Poisson 入力（各イベントで g += f_poi × w_syn）。

ギャップ結合は化学シナプスと別に、遅延 0 の「スパイクレット」（発火した瞬間に相手の v へ直接加算）で近似する。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numba import njit

from .data import Connectome


@dataclass
class LIFParams:
    dt: float = 0.1          # ms
    t_run: float = 1000.0    # ms
    v0: float = -52.0        # mV（静止電位）
    v_rst: float = -52.0
    v_th: float = -45.0
    tau_m: float = 20.0      # ms
    tau_syn: float = 5.0
    t_ref: float = 2.2
    t_delay: float = 1.8
    w_syn: float = 0.275     # mV / シナプス（Shiu 2024 が FlyWire で使った値）
    w_scale: float = 1.0     # データセット間の検出シナプス数の差を吸収する倍率
    # 受け手の大きさによる 1 シナプスの効きの補正:  × (size_ref / 受け手の総入力シナプス数) ** size_alpha
    #   alpha = 0 : Shiu 2024 と同じ（全細胞で 1 シナプスの効きが等しい）
    #   alpha = 1 : 入力割合で正規化（入力コンダクタンス ∝ 膜面積 ∝ 総シナプス数、という近似）
    # size_ref は MaleCNS の総入力シナプス数の中央値。中央値の細胞では alpha によらず重みが変わらない。
    size_alpha: float = 0.0
    size_ref: float = 346.0
    f_poi: float = 250.0     # Poisson 入力1発の重み（w_syn の倍数）
    adapt_b: float = 0.0     # mV。発火ごとの順応の増分
    tau_adapt: float = 100.0 # ms
    sign_preset: str = "fast_only"


@dataclass
class Network:
    """発火元(pre)ごとに並べた CSR。weights は mV。"""
    n: int
    indptr: np.ndarray
    indices: np.ndarray
    weights: np.ndarray
    gap_indptr: np.ndarray = field(default=None)
    gap_indices: np.ndarray = field(default=None)
    gap_weights: np.ndarray = field(default=None)

    @classmethod
    def from_connectome(cls, conn: Connectome, p: LIFParams, gaps=()):
        w = conn.signed_weights(p.sign_preset) * (p.w_syn * p.w_scale)
        if p.size_alpha:
            tot = np.bincount(conn.post, weights=conn.count, minlength=conn.n)
            w = (w * (p.size_ref / np.maximum(tot[conn.post], 1)) ** p.size_alpha).astype(np.float32)
        nz = w != 0
        net = cls(conn.n, *_csr(conn.pre[nz], conn.post[nz], w[nz], conn.n))
        return net.with_gaps(gaps)

    def with_gaps(self, gaps):
        """gaps: (pre_idx, post_idx, kick_mV) の並び。pre の発火で post の v に kick_mV を即時加算。"""
        if len(gaps):
            g = np.asarray(gaps, dtype=np.float64)
            ip, ix, gw = _csr(g[:, 0].astype(np.int32), g[:, 1].astype(np.int32), g[:, 2], self.n)
        else:
            ip, ix, gw = np.zeros(self.n + 1, np.int64), np.zeros(0, np.int32), np.zeros(0, np.float32)
        return Network(self.n, self.indptr, self.indices, self.weights, ip, ix, gw)


def _csr(pre, post, w, n):
    order = np.argsort(pre, kind="stable")
    indptr = np.zeros(n + 1, np.int64)
    np.add.at(indptr, pre + 1, 1)
    return np.cumsum(indptr), post[order].astype(np.int32), w[order].astype(np.float32)


@dataclass
class Result:
    counts: np.ndarray       # (trials, n) 発火数
    t_run: float             # ms
    spikes: list             # 試行ごとの (t_ms, idx) 配列。record 指定時のみ

    @property
    def rates(self) -> np.ndarray:
        """試行平均の発火率 [Hz]。"""
        return self.counts.mean(0) / (self.t_run / 1000.0)


def simulate(net: Network, stim=(), p: LIFParams = LIFParams(), n_trials: int = 10,
             silence=(), record=(), seed: int = 0) -> Result:
    """stim: [(idx 配列 or int, Poisson レート[Hz]), ...]。silence: 発火させない idx。record: スパイク時刻を残す idx。"""
    stim_idx, stim_rate = [], []
    for k, rate in stim:
        k = np.atleast_1d(k)
        stim_idx.append(k)
        stim_rate.append(np.full(len(k), rate, np.float64))
    stim_idx = np.concatenate(stim_idx).astype(np.int32) if stim_idx else np.zeros(0, np.int32)
    stim_p = np.concatenate(stim_rate) * p.dt / 1000.0 if stim_rate else np.zeros(0)
    silenced = np.zeros(net.n, np.bool_)
    silenced[np.asarray(silence, dtype=np.int64)] = True
    rec = np.zeros(net.n, np.bool_)
    rec[np.asarray(record, dtype=np.int64)] = True

    steps = int(round(p.t_run / p.dt))
    dly = max(1, int(round(p.t_delay / p.dt)))
    counts = np.zeros((n_trials, net.n), np.int32)
    spikes = []
    for tr in range(n_trials):
        c, st, si = _run(net.indptr, net.indices, net.weights,
                         net.gap_indptr, net.gap_indices, net.gap_weights,
                         net.n, steps, dly, stim_idx, stim_p, p.f_poi * p.w_syn,
                         silenced, rec, np.exp(-p.dt / p.tau_syn), p.dt / p.tau_m,
                         p.v0, p.v_th, p.v_rst, int(round(p.t_ref / p.dt)), seed * 100003 + tr,
                         p.adapt_b, np.exp(-p.dt / p.tau_adapt))
        counts[tr] = c
        if rec.any():
            spikes.append(np.stack([st * p.dt, si], 1))
    return Result(counts, p.t_run, spikes)


@njit(cache=True)
def _run(indptr, indices, weights, gptr, gidx, gw, n, steps, dly, stim_idx, stim_p, stim_w,
         silenced, rec, decay_g, dt_tau, v0, vth, vrst, ref_steps, seed, adapt_b, decay_a):
    np.random.seed(seed)
    v = np.full(n, v0)
    g = np.zeros(n)
    a = np.zeros(n)
    ref = np.zeros(n, np.int32)
    counts = np.zeros(n, np.int32)
    buf = np.empty((dly, n), np.int32)
    blen = np.zeros(dly, np.int64)
    cap = 1 << 16
    rec_t = np.empty(cap, np.int64)
    rec_i = np.empty(cap, np.int32)
    nrec = 0
    for step in range(steps):
        slot = step % dly
        # dly ステップ前の発火を配達
        for k in range(blen[slot]):
            pre = buf[slot, k]
            for j in range(indptr[pre], indptr[pre + 1]):
                g[indices[j]] += weights[j]
        blen[slot] = 0
        for s in range(stim_idx.shape[0]):
            if np.random.random() < stim_p[s]:
                g[stim_idx[s]] += stim_w
        for i in range(n):
            if adapt_b != 0.0:
                a[i] *= decay_a
            if ref[i] > 0:
                ref[i] -= 1
                continue
            v[i] += (v0 - v[i] + g[i] - a[i]) * dt_tau
            g[i] *= decay_g
        for i in range(n):
            if v[i] > vth and ref[i] == 0 and not silenced[i]:
                v[i] = vrst
                ref[i] = ref_steps
                a[i] += adapt_b
                counts[i] += 1
                buf[slot, blen[slot]] = i
                blen[slot] += 1
                for j in range(gptr[i], gptr[i + 1]):
                    v[gidx[j]] += gw[j]
                if rec[i]:
                    if nrec == cap:
                        cap *= 2
                        t2 = np.empty(cap, np.int64)
                        i2 = np.empty(cap, np.int32)
                        t2[:nrec] = rec_t[:nrec]
                        i2[:nrec] = rec_i[:nrec]
                        rec_t, rec_i = t2, i2
                    rec_t[nrec] = step
                    rec_i[nrec] = i
                    nrec += 1
    return counts, rec_t[:nrec], rec_i[:nrec]
