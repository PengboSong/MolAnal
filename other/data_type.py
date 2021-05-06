# coding=utf-8

import argparse

import numpy as np


GMXDataType = argparse.Namespace(
    BOOL=np.bool,
    INT=np.int32,
    UINT=np.uint32,
    LL=np.int64,
    ULL=np.uint64,
    FLOAT=np.float32,
    DOUBLE=np.float64,
    REAL=np.float64,
)
