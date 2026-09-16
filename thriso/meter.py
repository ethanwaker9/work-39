import time
from collections import defaultdict


class Meter:
    def __init__(self):
        self.reset()

    def reset(self):
        self.party = 0
        self.phase = "sign"
        self.ga = defaultdict(int)
        self.ga_party = defaultdict(lambda: defaultdict(int))
        self.round_ga = []
        self.curves = defaultdict(int)
        self.scalars = defaultdict(int)
        self.hashbytes = defaultdict(int)
        self.rawbytes = defaultdict(int)
        self.rounds = defaultdict(int)
        self.time = defaultdict(float)
        self._current_round = None
        self._t0 = None

    def set_party(self, i):
        self.party = i

    def set_phase(self, name):
        self.phase = name

    def count_ga(self, k=1):
        self.ga[self.phase] += k
        self.ga_party[self.phase][self.party] += k
        if self._current_round is not None:
            self._current_round[self.party] += k

    def begin_round(self):
        self._current_round = defaultdict(int)
        self.rounds[self.phase] += 1

    def end_round(self):
        if self._current_round is not None:
            self.round_ga.append((self.phase, dict(self._current_round)))
        self._current_round = None

    def send(self, curves=0, scalars=0, hashbytes=0, rawbytes=0):
        self.curves[self.phase] += curves
        self.scalars[self.phase] += scalars
        self.hashbytes[self.phase] += hashbytes
        self.rawbytes[self.phase] += rawbytes

    def start_clock(self):
        self._t0 = time.perf_counter()

    def stop_clock(self):
        if self._t0 is not None:
            self.time[self.phase] += time.perf_counter() - self._t0
        self._t0 = None

    def latency_ga(self, phase):
        total = 0
        for ph, rg in self.round_ga:
            if ph == phase and rg:
                total += max(rg.values())
        return total

    def bytes(self, phase, curve_bytes, scalar_bytes):
        return (self.curves[phase] * curve_bytes + self.scalars[phase] * scalar_bytes
                + self.hashbytes[phase] + self.rawbytes[phase])

    def summary(self, curve_bytes, scalar_bytes):
        phases = set(self.ga) | set(self.curves) | set(self.scalars) | set(self.hashbytes) | set(self.rawbytes) | set(self.rounds)
        out = {}
        for ph in sorted(phases):
            out[ph] = {
                "group_actions": self.ga[ph],
                "latency_group_actions": self.latency_ga(ph),
                "max_party_group_actions": max(self.ga_party[ph].values()) if self.ga_party[ph] else 0,
                "curves": self.curves[ph],
                "scalars": self.scalars[ph],
                "hash_bytes": self.hashbytes[ph],
                "raw_bytes": self.rawbytes[ph],
                "bytes": self.bytes(ph, curve_bytes, scalar_bytes),
                "rounds": self.rounds[ph],
                "seconds": self.time[ph],
            }
        return out


METER = Meter()
