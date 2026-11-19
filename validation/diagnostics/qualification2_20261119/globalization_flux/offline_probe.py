"""Offline access to the exact frozen residual. Invalid states raise explicitly."""
import ctypes as ct
from pathlib import Path
import os

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np

ROOT = Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119/offline_replay')
NX = 306
ND = np.ctypeslib.ndpointer(dtype=np.float64, flags='CONTIGUOUS')
FD = np.ctypeslib.ndpointer(dtype=np.float64, flags='F_CONTIGUOUS')
IP = ct.POINTER(ct.c_int)
DP = ct.POINTER(ct.c_double)


class InvalidState(RuntimeError):
    def __init__(self, status, index, stage):
        self.details = {'status': status, 'index': index, 'stage': stage}
        super().__init__(str(self.details))


class Probe:
    def __init__(self, case='D_SS1', mode=1):
        self.case = case
        self.mode = mode
        self.lib = ct.CDLL(str(ROOT/'build/offline_api.so'))
        self.lib.probe_init.argtypes = [ct.c_char_p, ct.c_int, ct.c_int, IP]
        self.lib.probe_state.argtypes = [ct.c_int, ND, ND, ND]
        self.lib.probe_eval.argtypes = [ND, ct.c_int, ND, FD, DP, ND, ND, IP, IP, IP]
        self.lib.probe_margins.argtypes = [ND, ND, IP]
        self.lib.probe_interface.argtypes = [ND, ND, IP]
        folder = str(ROOT/case).encode()
        code = ct.c_int()
        self.lib.probe_init(folder, len(folder), mode, ct.byref(code))
        if code.value:
            raise InvalidState(code.value, -1, 'PREPARE')
        self.old, self.scale, self.meta = self.state(0)
        self.trial, _, _ = self.state(1)
        self.calls = 0

    def state(self, which):
        x, scales, meta = np.empty(NX), np.empty(NX), np.empty(6)
        self.lib.probe_state(which, x, scales, meta)
        return x, scales, meta

    def evaluate(self, x, jac=True):
        self.calls += 1
        residual = np.empty(NX)
        matrix = np.empty((NX, NX), order='F')
        fin, fout = np.empty(3), np.empty(3)
        norm = ct.c_double()
        code, index, stage = ct.c_int(), ct.c_int(), ct.c_int()
        self.lib.probe_eval(np.ascontiguousarray(x), jac, residual, matrix, ct.byref(norm),
                            fin, fout, ct.byref(code), ct.byref(index), ct.byref(stage))
        if code.value:
            raise InvalidState(code.value, index.value, stage.value)
        return residual, matrix, norm.value, fin, fout

    def margins(self, x):
        margins = np.empty(53)
        code = ct.c_int()
        self.lib.probe_margins(np.ascontiguousarray(x), margins, ct.byref(code))
        if code.value:
            raise InvalidState(code.value, -1, 'MARGIN')
        return margins

    def interface(self, x):
        data = np.empty(48)
        code = ct.c_int()
        self.lib.probe_interface(np.ascontiguousarray(x), data, ct.byref(code))
        if code.value:
            raise InvalidState(code.value, -1, 'INTERFACE')
        return data
