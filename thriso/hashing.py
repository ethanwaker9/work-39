import hashlib
import secrets

LAMBDA_BYTES = 32


def _enc(x):
    if isinstance(x, bytes):
        return b"b" + len(x).to_bytes(4, "big") + x
    if isinstance(x, str):
        y = x.encode()
        return b"s" + len(y).to_bytes(4, "big") + y
    if isinstance(x, bool):
        return b"o" + (b"\x01" if x else b"\x00")
    if isinstance(x, int):
        sign = b"-" if x < 0 else b"+"
        y = abs(x)
        raw = y.to_bytes((y.bit_length() + 7) // 8 or 1, "big")
        return b"i" + sign + len(raw).to_bytes(4, "big") + raw
    if isinstance(x, (list, tuple)):
        body = b"".join(_enc(v) for v in x)
        return b"l" + len(x).to_bytes(4, "big") + body
    raise TypeError(type(x))


def xof(domain, *items, length=LAMBDA_BYTES):
    h = hashlib.shake_256()
    h.update(_enc(domain))
    for it in items:
        h.update(_enc(it))
    return h.digest(length)


def hash_to_int(domain, modulus, *items):
    raw = xof(domain, *items, length=(modulus.bit_length() + 7) // 8 + 16)
    return int.from_bytes(raw, "big") % modulus


def challenge_vector(domain, count, K, *items):
    size = 2 * K + 1
    out = []
    ctr = 0
    nbits = max(1, (size - 1).bit_length())
    nbytes = (nbits + 7) // 8 + 1
    while len(out) < count:
        raw = xof(domain, ctr, *items, length=nbytes * (count - len(out)) * 2)
        ctr += 1
        mask = (1 << nbits) - 1
        for j in range(0, len(raw), nbytes):
            v = int.from_bytes(raw[j:j + nbytes], "big") & mask
            if v < size:
                out.append(v - K)
                if len(out) == count:
                    break
    return out


def commit(value, rand=None):
    if rand is None:
        rand = secrets.token_bytes(LAMBDA_BYTES)
    return xof("commit", value, rand), rand


def prf(key, modulus, *items):
    return hash_to_int("prf", modulus, key, *items)


def random_key():
    return secrets.token_bytes(LAMBDA_BYTES)
