# coding=utf-8

from enum import Enum


class RegisterFunction(object):
    def __init__(self):
        self.register_table_ = []

    def register_(self, func_name, func, *args):
        """Register functions that lanuched with arguments from command-line"""
        self.register_table_[func_name] = (func, *args)

    def launch_(self, fullcmd):
        """Launch registered functions from command-line input"""
        fullcmd = fullcmd.strip()
        if ' ' in fullcmd:
            cmd, rawargs = fullcmd.split(' ', maxsplit=1)
        else:
            cmd, rawargs = fullcmd, ''
        regfunc = self.register_table_[cmd]
        # Split input argument string
        splitargs = []


        args = []
        for argtype in regfunc[1:]:
            pass