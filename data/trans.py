#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re

def main():
    re_name = re.compile('([OA])-(.+)@(.+)')
    
    f = open('basis2', 'r')
    for line in f:
        line = line.rstrip()

        matchObj = re_name.match(line)
        if matchObj:
            new_line = '%s-%s.%s' % (matchObj.group(1),
                                     matchObj.group(3),
                                     matchObj.group(2))
            #print('%s => %s' % (line, new_line))
            line = new_line

        print(line)
    f.close()


if __name__ == "__main__":
    main()
                
