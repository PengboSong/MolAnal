# coding=utf-8

from enum import Enum
from gmx.other.input_func import SplitArgs

class RegisterFunction(object):
    def __init__(self):
        self.register_table_ = []
        self.exit_signal_ = False

    def register_(self, func_name, func, *args):
        """Register functions that lanuched with arguments from command-line"""
        self.register_table_[func_name] = (func, *args)

    def launch_(self, fullcmd):
        """Launch registered functions from command-line input"""
        # Split input argument string
        splitargs = SplitArgs().digest(fullcmd)
        if len(splitargs) > 1:
            cmd = splitargs[0].lower()
            if cmd in self.register_table_:
                regfunc = self.register_table_[cmd][0]
                regfunc(*splitargs[1:])
            else:
                print("[Error] Unknown command {}. Available commands include:".format(cmd))
                for nm in self.register_table_:
                    print("- {}".format(nm))

    def help_(self):
        """Show help information"""
        for nm, detail in self.register_table_:
            print("* {}".format(nm))
            func = detail[0]
            print(func.__doc__)
    
    def exit_(self):
        """Close interactive console and exit"""
        self.exit_signal_ = True
