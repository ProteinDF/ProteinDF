#!/usr/bin/env python

import struct

class Endian:
    def __init__(self):
        num = struct.pack("@I", 0x01020304)
        self.big = (num== "\x01\x02\x03\x04")
        self.little = (num=="\x04\x03\x02\x01")

def main():
    e = Endian()
    if (e.big == True):
        print "This machine is big endian."
    elif (e.little == True):
        print "This machine is little endian."
    else:
        print "Sorry, endian of this machine is unknown."


if __name__ == '__main__':
    main()


