#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
from xml.parsers import expat

class Element(object):
    """parsed XML element"""
    def __init__(self, name, attributes = {}):
        self.name = name
        self.attributes = copy.deepcopy(attributes)
        
        self.cdata = ''
        self.children = []
    
    def add_child(self, element):
        self.children.append(element)

    def set_attribute(self, key, value):
        assert(isinstance(key, str))
        assert(isinstance(value, str))
        self.attributes[key] = value
        
    def get_attribute(self, key):
        return self.attributes.get(key, None)

    def get_data(self):
        return self.cdata

    def get_elements(self, name=''):
        if name:
            return [c for c in self.children if c.name == name]
        else:
            return list(self.children)

    def get_xml(self):
        xml = "<%s" % (self.name)
        for key in self.attributes.iterkeys():
            xml += ' %s="%s"' % (key, self.attributes[key])

        child_elements = self.get_elements()
        if (len(child_elements) != 0):
            xml += ">\n"
            for element in child_elements:
                xml += element.get_xml()
            xml += "</%s>\n" % (self.name)
        else:
            xml += " />\n"

        return xml;

    def __str__(self):
        return self.get_xml()


class Xml2Obj(object):
    """transpose Element from XML"""
    def __init__(self):
        self.root = None
        self.nodeStack = []

    def StartElement(self, name, attributes):
        """Expat handler of start element"""
        element = Element(name.encode(), attributes)
        if self.nodeStack:
            parent = self.nodeStack[-1]
            parent.add_child(element)
        else:
            self.root = element
        self.nodeStack.append(element)

    def EndElement(self, name):
        self.nodeStack.pop()

    def CharacterData(self, data):
        if data.strip():
            data = data.encode()
            element = self.nodeStack[-1]
            element.cdata += data
    
    def Parse(self, filename):
        Parser = expat.ParserCreate()
        Parser.StartElementHandler = self.StartElement
        Parser.EndElementHandler = self.EndElement
        Parser.CharacterDataHandler = self.CharacterData
        
        ParserStatus = Parser.Parse(open(filename).read(), 1)
        return self.root

def main():
    try:
        optlist, args = getopt.gnu_getopt(sys.argv[1:], "h", longopts=["help"])
    except getopt.GetoptError:
        # exit before display the help message
        usage()
        sys.exit(2)

    for opt, arg in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)

    if (len(args) == 0):
        usage()
        sys.exit(2)

    # setting
    xmlpath = args.pop(0)

    parser = Xml2Obj()
    root_element = parser.Parse(xmlpath)

    print root_element


if __name__ == '__main__':
    main()


