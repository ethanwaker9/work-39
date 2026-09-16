from .schemes.dm20 import DM20, DEMS24
from .schemes.sashimi import Sashimi, CSISharKActive
from .schemes.camposmuth import CamposMuth
from .schemes.thresher import ThreshERSharK
from .schemes.grass import GRASS, GRASSPlus
from .schemes.ours import Albacore

SCHEMES = {
    "dm20": DM20,
    "dems24": DEMS24,
    "sashimi": Sashimi,
    "csishark": CSISharKActive,
    "cm22": CamposMuth,
    "thresher": ThreshERSharK,
    "grass": GRASS,
    "grassplus": GRASSPlus,
    "ours": Albacore,
}

ORDER = ["dm20", "dems24", "sashimi", "csishark", "cm22", "thresher", "grass", "grassplus", "ours"]


def parties_for(key, t):
    if key in ("sashimi",):
        return t, t
    if key == "thresher":
        return 2 * t - 1, t
    if key == "grass":
        return t, t
    return 2 * t - 1, t


def signer_set(key, n, t):
    if key == "thresher":
        return list(range(1, n + 1))
    return list(range(1, t + 1))
