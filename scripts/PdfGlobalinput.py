#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re

class PdfGlobalinput(object):
    """
    parse and analize ProteinDF global input file (eg. fl_Globalinput)
    """

    # member variable ==================================================
    __file_path = ""
    __file = None # file object

    __current_group = ""
    __current_keyword = ""
    __current_value = ""
    __data = {}

    # comiled regulae explations =======================================
    __re_comment = re.compile("^(//|#)")
    __re_start_section = re.compile(">>>>\s*(\S+)\s*")
    __re_keyword_and_value = re.compile("\[(.*)\]\s+\[(.*)\]")

    #
    def __init__(self, file_path):
        if (os.access(file_path, os.R_OK) == True):
            self.__file_path = file_path
            self.__read()


    def __str__(self):
        output = ""
        for group in self.__data.keys():
            output += ">>>>%s\n" % (group)
            for keyword, value in self.__data[group].iteritems():
                output += "%s = %s\n" % (keyword, value)
        return output


    def get_groups(self):
        return self.__data.keys()


    def get_keywords(self, group):
        group = str(group)
        answer = None
        if (self.__data.has_key(group)):
            answer = self.__data[group].keys()
        return answer


    def get_value(self, group, keyword):
        group = str(group)
        keyword = str(keyword)
        answer = None
        if ((self.__data.has_key(group) == True) and
            (self.__data[group].has_key(keyword) == True)):
            answer = self.__data[group][keyword]
        return answer


    def __read(self):
        """read ProteinDF global input file (e.g. fl_Input/fl_Globalinput) and parse it."""

        self.__file = open(self.__file_path, "r")

        isInBracket = False
        while True:
            line = self.__file.readline()
            if (len(line) == 0):
                break
            line = line.rstrip('\n')

            if (self.__re_comment.search(line) != None):
                """"comment line"""
                continue
            
            while True:
                line = line.strip()
                #print line
                if (len(line) != 0):
                    if (self.__re_start_section.search(line) != None):
                        line = self.__read_start_section(line)
                        continue
                    elif (self.__re_keyword_and_value.search(line) != None):
                        line = self.__read_keyword_and_value(line)
                        continue
                    else:
                        print "parse error: " + line
                        sys.exit()
                else:
                    break

        self.__file.close()

    def __read_start_section(self, line):
        """set current section"""
        matchObj = self.__re_start_section.match(line)
        group = matchObj.group(1)
        self.__current_group = group
        return ""

    def __read_keyword_and_value(self, line):
        matchObj = self.__re_keyword_and_value.match(line)
        keyword = str(matchObj.group(1))
        value = str(matchObj.group(2))
        group = str(self.__current_group)
        if (self.__data.has_key(group) != True):
            self.__data[group] = {}
        self.__data[group][keyword] = value
        return ""


