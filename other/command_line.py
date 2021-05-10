# coding=utf-8

import re


class ConsoleBase(object):
    YES_PATTERN = re.compile(r"[yY][eE]?[sS]?")
    NO_PATTERN = re.compile(r"[nN][oO]?")

    def askyn(self, prompt):
        reply = input(prompt + "[Y(es)/N(o)]").strip()
        ymatch = re.match(self.YES_PATTERN, reply)
        nmatch = re.match(self.NO_PATTERN, reply)
        while not ymatch and not nmatch:
            reply = input(prompt + "[Y(es)/N(o)]").strip()
            ymatch = re.match(self.YES_PATTERN, reply)
            nmatch = re.match(self.NO_PATTERN, reply)
        if ymatch:
            return True
        else:
            return False
