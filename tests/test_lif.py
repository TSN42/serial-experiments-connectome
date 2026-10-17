import numpy as np

from secon.lif import LIFParams, Network, _csr, simulate


def chain(w_mv, n=3):
    pre = np.arange(n - 1, dtype=np.int32)
    post = pre + 1
    return Network(n, *_csr(pre, post, np.full(n - 1, w_mv), n)).with_gaps(())


P = LIFParams(t_run=200.0)


def test_unstimulated_network_is_silent():
    assert simulate(chain(20.0), [], P, n_trials=1).counts.sum() == 0


def test_strong_synapse_propagates_with_delay():
    r = simulate(chain(20.0), [(0, 150.0)], P, n_trials=1, record=[0, 1])
    t, i = r.spikes[0][:, 0], r.spikes[0][:, 1]
    assert r.counts[0, 0] > 10 and r.counts[0, 1] > 0
    assert t[i == 1].min() >= t[i == 0].min() + P.t_delay


def test_inhibitory_synapse_does_not_drive():
    assert simulate(chain(-20.0), [(0, 150.0)], P, n_trials=1).counts[0, 1] == 0


def test_silencing_blocks_propagation():
    r = simulate(chain(20.0), [(0, 150.0)], P, n_trials=1, silence=[1])
    assert r.counts[0, 1] == 0 and r.counts[0, 2] == 0


def test_gap_kick_is_immediate():
    net = chain(0.0).with_gaps([(0, 1, 10.0)])
    r = simulate(net, [(0, 150.0)], P, n_trials=1, record=[0, 1])
    t, i = r.spikes[0][:, 0], r.spikes[0][:, 1]
    assert r.counts[0, 1] > 0
    assert t[i == 1].min() - t[i == 0].min() < P.t_delay


def test_adaptation_lowers_rate():
    base = simulate(chain(20.0), [(0, 150.0)], P, n_trials=1).counts[0, 1]
    p = LIFParams(t_run=200.0, adapt_b=5.0)
    adapted = simulate(chain(20.0), [(0, 150.0)], p, n_trials=1).counts[0, 1]
    assert adapted < base
