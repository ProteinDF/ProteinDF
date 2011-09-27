#!/usr/bin/env python
# -*- coding: utf-8 -*-

import urllib
import urllib2
import re
import os
import optparse
import copy

#from HTMLParser import HTMLParser

#class MyHTMLParser(HTMLParser):
#    def handle_data(self, data):
#        data = data.strip(" \t\r\n")
#        if (len(data) > 0):
#            print 'data: "%s"' % (data)
#
#    def handle_starttag(self, tag, attrs):
#        #pass
#        print 'start tag: "%s"' % (tag)
#    
#    def handle_endtag(self, tag):
#        #pass
#        print 'end tag: "%s"' % (tag)

class BasisSetExchange(object):
    __list_url = 'https://bse.pnl.gov:443/bse/portal/user/anon/js_peid/11535052407933/panel/Main/template/content'
    __file_path = 'content.html'
    __basisset = {}
    __html = ""
    __data = {}

    def __init__(self, is_download = True):
        self.__basisset = {}
        self.__html = ""
        if (is_download == True):
            self.download(self.__list_url, self.__file_path)
            

    def download(self, url, file_path):
        print("downloading %s to %s... " % (url, file_path))
        proxies = {}
        opener = urllib.FancyURLopener(proxies)
        opener.retrieve(url, file_path)
        print("finished.\n")


    def wget(self, url, file_path):
        cmd = 'wget -O %s "%s"' % (file_path, url)
        os.system(cmd)


    def read_text(self, file_path):
        print("reading %s..." % (file_path))
        fr = open(file_path, 'r')
        html = fr.read()
        fr.close()
        print("finished.\n")

        return html


    def listup_basisset(self, isPrinting = False):
        html = self.read_text(self.__file_path)
        print("listing...\n")
        re_basisset = re.compile('basisSets\[(\d+)\] = new basisSet\("(.*?)",\s*"(.*?)",\s*"(.*?)",\s*"\[(.*?)\]",\s*"(.*?)",\s*"(.*?)",.*\)',
                                 re.MULTILINE)
        re_elements = re.compile('[,\s]+')
        for match_obj in re_basisset.finditer(html):
            index = int(match_obj.group(1))
            url = match_obj.group(2)
            name = match_obj.group(3)
            type = match_obj.group(4)
            elements = match_obj.group(5)
            status = match_obj.group(6)
            hasECP = match_obj.group(7)
            self.__basisset.setdefault(index, {})
            item = {}
            item['url'] = url
            item['name'] = name
            item['type'] = type
            item['elements'] = re_elements.split(elements)
            item['status'] = status
            item['hasECP'] = hasECP
            self.__basisset[index] = item
            if (isPrinting == True):
                print '[%d] "%s"' % (index, name)

    def download_basisset(self, index):
        elements = ""
        for e in self.__basisset[index]['elements']:
            elements += e + ' '
        elements = elements.strip()

        params = {}
        params['bsurl'] = self.__basisset[index]['url']
        params['bsname'] = self.__basisset[index]['name']
        params['elts'] = elements
        params['format'] = 'Gaussian94'
        params['minimize'] = 'checked'
        
        url = "https://bse.pnl.gov:443"
        url += "/bse/portal/user/anon/js_peid/11535052407933/action/portlets.BasisSetAction"
        url += "/template/courier_content/panel/Main"
        url += "/eventSubmit_doDownload/true?"
        url += urllib.urlencode(params)

        file_path = "basisset_%d.html" % (index)
        #print url
        #self.download(url, file_path)
        if (os.path.exists(file_path) == False):
            print("downloading...: %s" % (url))
            self.wget(url, file_path)
        else:
            print("file is existed. no downloaging.: %s" % (file_path))

        
    def parse_basisset_html(self, index):
        index = int(index)
        file_path = "basisset_%d.html" % (index)

        re_atom = re.compile('(\S\S?)\s+0')
        re_shell = re.compile('(\S+)\s+(\d+)\s+(\S+)')
        re_param = re.compile('\s*([+-]?\d\S*)\s+([+-]?\d\S+)')

        name = self.__basisset[index]['name']
        self.__data.setdefault(name, {})

        state = ''
        atom = ''
        orb_type = ''
        items = 0
        items_counter = 0
        shell_type = ''
        cGTO = []
        atomEncount = {} # key is atom
        f = open(file_path, 'r')
        for line in f:
            line = line.rstrip()
            #print line
            
            if (state == 'finished'):
                break

            if (line == '****'):
                atom = ''
                orb_type = ''
                items = 0
                items_counter = 0
                shell_type = ''
                cGTO = []
                state = 'begin'
                continue

            if (state == ''):
                continue

            if (state == 'begin'):
                match = re_atom.match(line)
                if (match != None):
                    atom = match.group(1)
                    atom = atom.lower()
                    atom = atom.capitalize()
                    if (atomEncount.has_key(atom) == False):
                        atomEncount[atom] = 0
                        self.__data[name].setdefault(atom, {})
                    else:
                        atomEncount[atom] = atomEncount[atom] +1;
                    if (atomEncount[atom] == 0):
                        orb_type = "O"
                    else:
                        orb_type = "A%d" % (atomEncount[atom])
                    self.__data[name][atom].setdefault(orb_type, [])
                    state = 'shell'
                    continue
                else:
                    state = ''
                    continue
            elif (state == 'shell'):
                if (len(atom) != 0):
                    match = re_shell.match(line)
                    assert(match != None)
                    shell_type = match.group(1)
                    items = int(match.group(2))
                    state = 'function'
                    continue
            elif (state == 'function'):
                match = re_param.match(line)
                assert(match != None)
                pGTO = {}
                pGTO['exponent'] = float(match.group(1))
                pGTO['coef'] = float(match.group(2))
                cGTO.append(pGTO)
                if (len(cGTO) == items):
                    state = 'end'

            if (state == 'end'):
                cgto_set = {}
                cgto_set['shell_type'] = shell_type
                cgto_set['cGTO'] = cGTO
                cGTO = []
                self.__data[name][atom][orb_type].append(cgto_set)
                state = 'shell'

        f.close()


    def __str__(self):
        output = ''
        for name in self.__data.keys():
            for atom in self.__data[name].keys():
                for orb_type in self.__data[name][atom].keys():
                    output += '%s-%s@%s\n' % (orb_type, atom, name)
                    for cgto_set in (self.__data[name][atom][orb_type]):
                        shell_type = cgto_set['shell_type']
                        cGTO = cgto_set['cGTO']
                        output += 'shell = %s\n' % (shell_type)
                        for pGTO in (cGTO):
                            exponent = pGTO['exponent']
                            coef = pGTO['coef']
                            output += '  %f %f\n' % (exponent, coef)
                    output += '\n\n'

        return output


    def output_pdf_format(self):
        output = ''
        for name in self.__data.keys():
            for atom in self.__data[name].keys():
                basis_name = '%s@%s' % (atom, name)

                if (self.__data[name][atom].has_key("O") == True):
                    output += 'O-%s\n' % (basis_name)
                    output += self.output_pdf_format2(name, atom, "O")
                    output += "\n"
                if ((self.__data[name][atom].has_key("A1") == True) and
                    (self.__data[name][atom].has_key("A2") == True)):
                    output += 'A-%s\n' % (basis_name)
                    output += self.output_pdf_format2(name, atom, "A1")
                    output += self.output_pdf_format2(name, atom, "A2")
                    output += "\n"
        return output

        
    def output_pdf_format2(self, name, atom, orb_type):
        shell_type_list = ('S', 'P', 'D')
        output = ""

        # expand shell type: "SP" ...
        max_index = len(self.__data[name][atom][orb_type])
        for index in range(max_index):
            cgto_set = self.__data[name][atom][orb_type][index]
            shell_type = cgto_set['shell_type']
            if (shell_type == "SP"):
                cgto_set2 = copy.deepcopy(cgto_set)
                cgto_set['shell_type'] = "S"
                cgto_set2['shell_type'] = "P"
                self.__data[name][atom][orb_type][index] = cgto_set
                self.__data[name][atom][orb_type].append(cgto_set2)

        # count up shell type
        s_type = 0
        p_type = 0
        d_type = 0
        f_type = 0
        g_type = 0
        i_type = 0
        for cgto_set in (self.__data[name][atom][orb_type]):
            shell_type = cgto_set['shell_type']
            if (shell_type == 'S'):
                s_type += 1
            elif (shell_type == 'P'):
                p_type += 1
            elif (shell_type == 'D'):
                d_type += 1
            elif (shell_type == 'F'):
                f_type += 1
            elif (shell_type == 'G'):
                g_type += 1
            elif (shell_type == 'I'):
                i_type += 1
            else:
                print("unknown shell_type: %s" % (shell_type))

        #if ((f_type > 0) or (g_type > 0) or (i_type > 0)):
        #    print 'sorry unsupported: %s\n' % (basis_name)
        #    continue

        if (orb_type == "O"):
            output += '    %d %d %d\n' % (s_type, p_type, d_type)
        else:
            output += '    %d %d %d 0\n' % (s_type, p_type, d_type)

        for cgto_set in (self.__data[name][atom][orb_type]):
            shell_type = cgto_set['shell_type']
            for shell in shell_type_list:
                if (shell_type != shell):
                    continue

                cGTO = cgto_set['cGTO']
                output += '    %d\n' % (len(cGTO))

                for pGTO in (cGTO):
                    exponent = pGTO['exponent']
                    coef = pGTO['coef']
                    if (orb_type == "O"):
                        output += '    %16.8f\t%16.8f\n' % (exponent, coef)
                    else:
                        output += '    %16.8f\n' % (exponent)
                #output += '\n'
        return output

        
def main():
    # parse options
    #
    parser = optparse.OptionParser(usage = "%prog [options]", version = "%prog 1.0")
    # set option
    parser.add_option("-d", "--download", dest="isDownload",
                      help="download basis set list",
                      action="store_true", default=False)

    parser.add_option("-w", "--write", dest="mpac_path",
                      help="output MsgPack file (default: results.mpac)", metavar="FILE",
                      default="results.mpac")
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="verbose",
                      action="store_true", default=False)
    (opts, args) = parser.parse_args()

    # setup
    isDownload = opts.isDownload
    isListup = False
    basisSetIndex = -1
    if (len(args) == 0):
        isListup = True
    else:
        basisSetIndex = int(args.pop())

    # main =====================================================================
    bse = BasisSetExchange(isDownload)
    bse.listup_basisset(isListup)

    if (basisSetIndex >= 0):
        bse.download_basisset(basisSetIndex)
        bse.parse_basisset_html(basisSetIndex)
        print bse.output_pdf_format()

if __name__ == '__main__':
    main()

