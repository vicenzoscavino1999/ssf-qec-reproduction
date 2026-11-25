#!/usr/bin/env python3
"""T2-2, piloto m = 6 (BORRADOR_PLAN_T2_2.txt v3, 2.11): registro de bloques, registros JSONL y manifiesto.

Codigo preparado para revision (version 1). Sin dependencias del proyecto.

Formato binario de un archivo de cadena (uno por nivel, semilla y repeticion), solo anadiendo:
- cabecera de bloque, 32 bytes, little-endian ('<4sBBBBIHHIIII'):
    magia b'T22B', version 1, nivel, semilla, repeticion (u8 cada uno), indice de bloque (u32),
    entradas (u16), muestras (u16), aceptaciones por peso, llamadas, llamadas evitadas, transiciones (u32);
  los pasos del bloque son sus muestras (una muestra por paso, sin E_0);
- entradas del histograma conjunto (firma, multiplicidad), en orden canonico de firma:
    multiplicidad (u32), d (u8) y d pares (clave u16, recuento u8): 5 + 3d bytes.
Validacion de rangos ANTES de escribir (C2): un campo fuera de rango lanza FormatoError y no se escribe nada.
Ademas: claves estrictamente crecientes y < numero de claves, recuentos >= 1, suma de multiplicidades =
muestras, aceptaciones = llamadas + evitadas, transiciones <= llamadas. Que m = 6 quepa no autoriza m = 48.
El lector valida lo mismo, la secuencia de bloques (0, 1, 2, ...) y que no sobren ni falten bytes.
"""
import hashlib
import json
import os
import struct

MAGIC = b"T22B"
VERSION = 1
HDR = struct.Struct("<4sBBBBIHHIIII")
assert HDR.size == 32
ENT = struct.Struct("<IB")
PAIR = struct.Struct("<HB")
U8, U16, U32 = 0xFF, 0xFFFF, 0xFFFFFFFF


class FormatoError(ValueError):
    """Campo fuera del rango del formato o archivo mal formado."""


class DatosNoCaben(Exception):
    """DATOS_NO_CABEN (7.1): lo escrito alcanzaria el techo de datos."""


def _rng(name, v, hi, lo=0):
    if not isinstance(v, int) or isinstance(v, bool) or not lo <= v <= hi:
        raise FormatoError("%s = %r fuera de rango [%d, %d]" % (name, v, lo, hi))


def encode_block(meta, hist, nkeys):
    """meta: nivel, semilla, repeticion, bloque, muestras, aceptaciones, llamadas, evitadas, transiciones.
    hist: dict firma -> multiplicidad; firma = tupla de (clave, recuento) con claves crecientes."""
    _rng("nivel", meta["nivel"], U8)
    _rng("semilla", meta["semilla"], U8)
    _rng("repeticion", meta["repeticion"], U8)
    _rng("bloque", meta["bloque"], U32)
    _rng("muestras", meta["muestras"], U16, 1)
    for k in ("aceptaciones", "llamadas", "evitadas", "transiciones"):
        _rng(k, meta[k], meta["muestras"])
    if meta["aceptaciones"] != meta["llamadas"] + meta["evitadas"]:
        raise FormatoError("aceptaciones != llamadas + evitadas")
    if meta["transiciones"] > meta["llamadas"]:
        raise FormatoError("transiciones > llamadas")
    _rng("numero de claves", nkeys, U16 + 1, 1)
    _rng("entradas", len(hist), U16, 1)
    parts = []
    tot = 0
    for sig in sorted(hist):
        m = hist[sig]
        _rng("multiplicidad", m, U32, 1)
        _rng("d", len(sig), U8)
        prev = -1
        body = [ENT.pack(m, len(sig))]
        for k, n in sig:
            _rng("clave", k, nkeys - 1)
            if k <= prev:
                raise FormatoError("claves no crecientes en una firma")
            prev = k
            _rng("recuento", n, U8, 1)
            body.append(PAIR.pack(k, n))
        parts.append(b"".join(body))
        tot += m
    if tot != meta["muestras"]:
        raise FormatoError("suma de multiplicidades %d != muestras %d" % (tot, meta["muestras"]))
    hdr = HDR.pack(MAGIC, VERSION, meta["nivel"], meta["semilla"], meta["repeticion"], meta["bloque"],
                   len(hist), meta["muestras"], meta["aceptaciones"], meta["llamadas"], meta["evitadas"],
                   meta["transiciones"])
    return hdr + b"".join(parts)


def decode_file(raw, nkeys, expect=None):
    """Lee todos los bloques de un archivo de cadena. expect: (nivel, semilla, repeticion) o None."""
    off, blocks = 0, []
    n = len(raw)
    while off < n:
        if n - off < HDR.size:
            raise FormatoError("cabecera truncada en el byte %d" % off)
        (mg, ver, lv, sd, rp, bi, ne, ns, acc, calls, av, tr) = HDR.unpack_from(raw, off)
        off += HDR.size
        if mg != MAGIC or ver != VERSION:
            raise FormatoError("magia o version distintas en el bloque %d" % len(blocks))
        if bi != len(blocks):
            raise FormatoError("bloque %d fuera de secuencia (se esperaba %d)" % (bi, len(blocks)))
        if expect is not None and (lv, sd, rp) != tuple(expect):
            raise FormatoError("bloque %d de otra cadena: %r" % (bi, (lv, sd, rp)))
        if ne < 1 or ns < 1 or acc != calls + av or tr > calls or acc > ns:
            raise FormatoError("contadores incoherentes en el bloque %d" % bi)
        hist, tot, prevsig = {}, 0, None
        for _ in range(ne):
            if n - off < ENT.size:
                raise FormatoError("entrada truncada en el bloque %d" % bi)
            m, d = ENT.unpack_from(raw, off)
            off += ENT.size
            if m < 1:
                raise FormatoError("multiplicidad nula en el bloque %d" % bi)
            if n - off < d * PAIR.size:
                raise FormatoError("pares truncados en el bloque %d" % bi)
            sig, prev = [], -1
            for _ in range(d):
                k, c = PAIR.unpack_from(raw, off)
                off += PAIR.size
                if k <= prev or k >= nkeys or c < 1:
                    raise FormatoError("par clave-recuento no valido en el bloque %d" % bi)
                prev = k
                sig.append((k, c))
            sig = tuple(sig)
            if prevsig is not None and not sig > prevsig:
                raise FormatoError("firmas fuera del orden canonico en el bloque %d" % bi)
            prevsig = sig
            hist[sig] = m
            tot += m
        if tot != ns:
            raise FormatoError("multiplicidades que no suman las muestras en el bloque %d" % bi)
        blocks.append(({"nivel": lv, "semilla": sd, "repeticion": rp, "bloque": bi, "muestras": ns,
                        "aceptaciones": acc, "llamadas": calls, "evitadas": av, "transiciones": tr}, hist))
    return blocks


class Budget:
    """Reserva de bytes en un contador compartido (multiprocessing.Value 'q') con techo estricto."""

    def __init__(self, counter, cap):
        self.c, self.cap = counter, cap

    def reserve(self, nbytes):
        with self.c.get_lock():
            if self.c.value + nbytes > self.cap:
                raise DatosNoCaben("el registro alcanzaria el techo de datos (%d + %d > %d bytes)"
                                   % (self.c.value, nbytes, self.cap))
            self.c.value += nbytes


class BlockWriter:
    """Escritura solo anadiendo de los bloques de una cadena, con SHA-256 acumulado. El archivo no debe
    existir. Cada bloque se codifica y valida entero antes de reservar bytes y escribir."""

    def __init__(self, path, key, nkeys, budget=None):
        if os.path.exists(path):
            raise FormatoError("el archivo de cadena ya existe: %s" % path)
        self.path, self.key, self.nkeys, self.budget = path, key, nkeys, budget
        self.h = hashlib.sha256()
        self.bytes = 0
        self.nblocks = 0
        self.fh = open(path, "ab")

    def write(self, counters, hist):
        meta = dict(counters, nivel=self.key[0], semilla=self.key[1], repeticion=self.key[2],
                    bloque=self.nblocks)
        raw = encode_block(meta, hist, self.nkeys)
        if self.budget is not None:
            self.budget.reserve(len(raw))
        self.fh.write(raw)
        self.fh.flush()
        self.h.update(raw)
        self.bytes += len(raw)
        self.nblocks += 1
        return len(raw)

    def close(self):
        if self.fh is not None:
            self.fh.flush()
            os.fsync(self.fh.fileno())
            self.fh.close()
            self.fh = None
        return {"archivo": os.path.basename(self.path), "sha256": self.h.hexdigest(), "bytes": self.bytes,
                "bloques": self.nblocks}


def append_jsonl(path, obj):
    line = json.dumps(obj, ensure_ascii=True, sort_keys=True, separators=(",", ":")) + "\n"
    with open(path, "a", encoding="ascii", newline="\n") as fo:
        fo.write(line)
        fo.flush()
        os.fsync(fo.fileno())
    return len(line)


def read_jsonl(path):
    out = []
    with open(path, encoding="ascii") as fh:
        for i, line in enumerate(fh):
            if not line.endswith("\n"):
                raise FormatoError("%s: linea %d sin terminar" % (path, i + 1))
            out.append(json.loads(line))
    return out


def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def write_json(path, obj, overwrite=False):
    """JSON ASCII con LF; escritura atomica (archivo temporal y os.replace). No sobrescribe salvo orden."""
    if os.path.exists(path) and not overwrite:
        raise FileExistsError(path)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="ascii", newline="\n") as fo:
        fo.write(json.dumps(obj, ensure_ascii=True, indent=1, sort_keys=True) + "\n")
        fo.flush()
        os.fsync(fo.fileno())
    os.replace(tmp, path)
    return sha_file(path)


def manifest(root, exclude=()):
    """Lista ordenada (ruta relativa con '/', sha256, bytes) de todos los archivos bajo root."""
    out = []
    for d, _, fs in os.walk(root):
        for f in fs:
            p = os.path.join(d, f)
            rel = os.path.relpath(p, root).replace("\\", "/")
            if rel in exclude:
                continue
            out.append([rel, sha_file(p), os.path.getsize(p)])
    return sorted(out)


def dir_bytes(root):
    tot = 0
    for d, _, fs in os.walk(root):
        for f in fs:
            try:
                tot += os.path.getsize(os.path.join(d, f))
            except OSError as ex:
                raise OSError("no se pudo medir %s: %s" % (f, ex))
    return tot
