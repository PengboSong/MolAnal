# coding = utf-8

import os

def valid_filetype(fpath, filetype="pdb"):
    filetype = '.' + filetype.lower()
    return os.path.isfile(fpath) and os.path.splitext(fpath)[-1] == filetype

def valid_onlytype(fpath, filetype="pdb"):
    filetype = '.' + filetype.lower()
    return os.path.splitext(fpath)[-1] == filetype
