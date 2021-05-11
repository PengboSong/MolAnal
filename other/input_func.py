# coding=utf-8

import collections
import re


def is_int(string):
    """Check whether a string can be converted to an integer"""
    assert isinstance(
        string, str), "Function is_int can only be used to check string."
    if re.match(r"^[\+\-]? *[0-9]+$", string.strip()):
        return True
    else:
        return False


def is_float(string):
    """Check whether a string can be converted to a floating number"""
    assert isinstance(
        string, str), "Function is_float can only be used to check string."
    floatp1 = re.compile(r"^[\+\-]? *[0-9]+\.[0-9]*[eE]?[\+\-]?[0-9]*$")
    floatp2 = re.compile(r"^[\+\-]? *[0-9]*\.[0-9]+[eE]?[\+\-]?[0-9]*$")
    string = string.strip()
    if re.match(floatp1, string) or re.match(floatp2, string):
        return True
    else:
        return False


class SplitArgs(object):
    BRACKETS = {'(': ')', '[': ']', '{': '}'}
    QUOTES = ["'", '"']

    def __init__(self):
        self.quote_stack = []
        self.quote_buf = []
        self.stringbuf = []

    def clear(self):
        self.quote_stack.clear()
        self.quote_buf.clear()
        self.stringbuf.clear()

    def convert(self, string):
        """Convert string to integer, floating number if possible"""
        if is_int(string):
            return int(string)
        elif is_float(string):
            return float(string)
        else:
            return string.strip()

    def convert_list(self, inlist, closed_char):
        """Convert list to correct Python object determined by closed character"""
        if closed_char == ')':
            return tuple(inlist)
        elif closed_char == ']':
            return list(inlist)
        elif closed_char == '}':
            if all([isinstance(item, str) and ':' in item for item in inlist]):
                outdict = {}
                for item in inlist:
                    k, v = item.split(':', maxsplit=1)
                    outdict[self.convert(k)] = self.convert(v)
                return outdict
            else:
                return set(inlist)

    def digest(self, string):
        """Digest input string and convert it to Python object"""
        # Reverse string characters order to pop character from last
        self.stringbuf = list(string.strip()[::-1])
        return self.handle()

    def handle(self):
        """Convert substring to Python objects"""
        # Put splitted arguments into list
        obj = []
        # String buffer for each arguments
        argbuf = ''
        # Single character
        char = ''
        # Separator, set to space if at root level, otherwise set to semicolon
        sep = ' ' if len(self.quote_stack) == 0 else ','
        # Append converted Python object to container
        def objapp(obj, arg):
            if len(arg) != 0:
                obj.append(self.convert(arg))
                arg = ''
            return arg

        # Loop until string is digested
        while len(self.stringbuf) != 0:
            char = self.stringbuf.pop()
            if char == sep:   # Split arguments when hitting separator
                argbuf = objapp(obj, argbuf)
            elif char in self.BRACKETS:   # Hit left bracket
                self.quote_stack.append(self.BRACKETS[char])
                handle_obj = self.handle()
                if handle_obj:
                    obj.append(handle_obj)
            elif char in self.BRACKETS.values():   # Hit right bracket
                if len(self.quote_stack) == 0:   # Missing right bracket
                    raise ValueError(
                        "Found non-closed paired character. Invalid input string.")
                elif char != self.quote_stack[-1]:   # Bracket scope overlaps
                    raise ValueError(
                        "Found dismatched paired character. Invalid input string.")
                else:   # Paired right bracket
                    argbuf = objapp(obj, argbuf)
                    return self.convert_list(obj, self.quote_stack.pop())
            elif char in self.QUOTES:   # Hit quotes
                self.quote_stack.append(char)
                # Put characters into container until hitting another paired quote
                char = self.stringbuf.pop()
                while len(self.stringbuf) != 0 and char != self.quote_stack[-1]:
                    argbuf += char
                    char = self.stringbuf.pop()
                self.quote_stack.pop()
                argbuf = objapp(obj, argbuf)
            else:   # Hit regular character
                argbuf += char
        argbuf = objapp(obj, argbuf)
        if len(self.quote_stack) != 0:
            raise ValueError(
                "Found non-closed paired character. Invalid input string.")
        return obj


def convert_input_type(string):
    """Convert string to Python object(s).

    Splitted results are stored in a list. If list is empty, returns None.
    If list has only one element, returns this element. Otherwise, returns
    a full list with arguments in order.
    """
    res = SplitArgs().digest(string)
    if len(res) == 0:
        return None
    elif len(res) == 1:
        return res[0]
    else:
        return res


def verify_selection(ls, desc):
    if isinstance(ls, (list, tuple, set)):
        for v in ls:
            print('- {}'.format(v))
    elif isinstance(ls, dict):
        for k, v in ls.items():
            print('{} : {}'.format(k, v))
    else:
        raise TypeError("Invalid selection list."
                        "Valid types include list, set, tuple and dict.")
    
    inp = convert_input_type(input(desc))
    while inp not in ls:
        print("[Info] Got invalid option. Type again.")
        inp = convert_input_type(input(desc))
    return inp