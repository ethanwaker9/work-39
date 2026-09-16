# Actively Secure $`(t,n)`$-Threshold Isogeny Signatures

This repository contains our research implementation of *Albacore: Efficient Actively Secure Threshold Isogeny Signatures*. It implements

* **lifted Shamir sharing** over the class group order `N`, which supports any number of parties without restricting the key space,
* **Albacore**, a `(t, n)` threshold CSI-FiSh signature with committee key chains and audited signing chains, which needs no zero-knowledge proof during signing, and
* **eight previous threshold signatures from isogeny group actions**, implemented from their specifications in the same framework, together with the benchmark, figure and table scripts used in our research.

## Repository Layout

```
final_code/
├── thriso/                    core library
│   ├── csidh.py               CSIDH group action with Velu isogenies in twisted Edwards form
│   ├── classgroup.py          binary quadratic forms, class number, discrete logarithms, LLL, Babai
│   ├── groupaction.py         real toy CSIDH backend, counting backend, CSIDH-512 timing backend
│   ├── galois.py              Galois rings GR(p^e, d)
│   ├── sharing.py             Shamir over subgroups, lifted Shamir, additive, replicated, recursive sharing
│   ├── nizk.py                sigma protocols (Fiat-Shamir and Unruh), piecewise verifiable proofs, VSS proofs
│   ├── hashing.py             SHAKE256 based oracles, commitments and challenge derivation
│   ├── meter.py               counting of group actions, rounds, latency and bytes
│   ├── registry.py            scheme registry and party configurations
│   ├── toy_params.json        cached class group data of the toy CSIDH prime
│   └── schemes/
│       ├── base.py            CSI-FiSh with K public curves, parameter sets
│       ├── dm20.py            DM20 and DEMS24
│       ├── sashimi.py         CS20 (Sashimi) and ABCP23a (actively secure CSI-SharK signing)
│       ├── camposmuth.py      CM22
│       ├── thresher.py        ABCP23b (ThreshER SharK)
│       ├── grass.py           BBMP24 (GRASS) and BBDMP25 (GRASS+)
│       └── ours.py            Albacore
├── bench/
│   ├── run_bench.py           all experiments of our research
│   ├── make_tables.py         formula validation 
│   └── make_figures.py        EPS results
├── tests/
│   └── test_all.py            unit and integration tests
├── results/                   outputs of the experiments used in our research
```

## Implemented schemes

| Key         | Label    | Scheme                                                    | Sharing                         | Security          |
|-------------|----------|-----------------------------------------------------------|---------------------------------|-------------------|
| `dm20`      | DM20     | Threshold CSI-FiSh of De Feo and Meyer                    | Shamir over a subgroup          | passive           |
| `dems24`    | DEMS24   | Rerandomizable threshold CSI-FiSh of Das et al.           | Shamir over a subgroup          | passive           |
| `sashimi`   | CS20     | Sashimi of Cozzo and Smart                                 | additive, `n = t`               | active with abort |
| `csishark`  | ABCP23a  | Actively secure signing with CSI-SharK keys               | Shamir, structured keys         | active with abort |
| `cm22`      | CM22     | Campos and Muth                                            | Shamir with piecewise proofs    | robust            |
| `thresher`  | ABCP23b  | ThreshER SharK of Atapoor et al.                           | Shamir, honest majority         | robust            |
| `grass`     | BBMP24   | GRASS of Battagliola et al.                                | additive or replicated          | active with abort |
| `grassplus` | BBDMP25  | GRASS+ of Battagliola et al.                               | recursive replicated            | active, adaptive  |
| `ours`      | Albacore | This work                                                  | lifted Shamir                   | active with abort |

## Group action backends

* **Real toy CSIDH backend** (`ToyCSIDHAction`): the CSIDH action over the 56 bit prime `p = 4 * 3 * 5 * ... * 47 - 1` with class number `N = 405294675 = 3 * 5^2 * 61 * 88589`. Scalars are reduced with the lattice of relations and Babai's nearest plane algorithm as in CSI-FiSh. Since 3 and 5 divide `N`, the lifted sharing uses nontrivial Galois ring extensions on this backend.
* **Counting backend** (`MockAction`): represents a curve by its exponent modulo the CSIDH-512 class number and counts each evaluation of the action. It runs the unchanged protocol logic, including all proofs and checks, with CSIDH-512 sizes of 64 bytes per curve and 33 bytes per scalar.
* **CSIDH-512 timing backend** (`CSIDH512Timer`): evaluates the CSIDH-512 action on exponent vectors in `[-5, 5]^74` and measures the time of one group action.


## Running the Experiments

```bash
python3 tests/test_all.py -v
```
The tests check the compatibility of the real CSIDH action, the correctness of lifted Shamir sharing for many parties and with the point at infinity, the exact privacy of lifted sharing on small parameters, the validity of signatures of all nine schemes on the counting backend and of DM20 and Albacore on the real isogeny backend, the agreement of the counted group actions with the formulas of our research, and the detection of deviations by the audits of Albacore. They take about 25 seconds.

```bash
python3 bench/run_bench.py
python3 bench/make_tables.py
python3 bench/make_figures.py
```

`run_bench.py` writes `results/bench.json` after each stage and prints progress to the terminal. The full run takes about 50 minutes on one core, most of it for the nine schemes on the real isogeny backend. `python3 bench/run_bench.py --quick` runs a reduced configuration in a few minutes, and `--skip-real` omits the real backend.

The benchmark consists of
1. the time of one CSIDH-512 group action (20 samples),
2. lifted and subgroup share sizes for CSIDH-512 up to `n = 10^9`,
3. signing costs of all schemes for `t in {2, 3, 4, 5, 6, 8, 10, 12, 16}` with `K = 16`,
4. signing costs for `K in {1, 16, 256}` with `T in {71, 23, 13}`,
5. signing costs for proof repetitions `(R_G, R_S) in {(112, 71), (256, 162)}`,
6. GRASS with replicated sharing for `n = 2t - 1`,
7. 600 honest sessions and 600 sessions per deviation type of Albacore with `t = 4`, `n = 7`,
8. all schemes on the real CSIDH backend for `t in {2, 3}`.

`make_tables.py` verifies that each counted number of group actions equals the formulas of the complexity analysis and exits with an error otherwise.

## Library Usage

```python
from thriso.groupaction import MockAction, ToyCSIDHAction
from thriso.schemes.ours import Albacore

ga = ToyCSIDHAction()
scheme = Albacore(ga, n=5, t=3, K=16)
scheme.keygen()
committee = [2, 4, 5]
scheme.setup(committee)
signature = scheme.sign(b"message", committee)
assert scheme.verify(b"message", signature)
```

Every scheme class offers `keygen()`, `setup(S)`, `sign(message, S)` and `verify(message, signature)`, and the global meter `thriso.meter.METER` gives the number of group actions, the latency, the rounds and the bytes of each phase.

## Notes
The implementation is intended for measurement and validation. It is not constant time and must not be used to protect real keys. All parties of a protocol run in one process, and the meter attributes each group action and each message to the party that computes or sends it. The experiments measure committee setup and signing. Key generation runs as a trusted dealer for all schemes except BBMP24, whose key generation chain with proofs is executed as specified; for ABCP23b and BBDMP25 the dealer produces the same share structure as their distributed or recursive key generation.
