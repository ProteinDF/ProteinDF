#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import optparse
import re
#import zipfile
import tarfile

from PdfGlobalinput import PdfGlobalinput

PDF_FL_INPUT = "fl_Input"
PDF_FL_WORK = "fl_Work"
PDF_GLOBALINPUT = "fl_Globalinput"

alive_list = ('work-up.zip',
              'QcStep',
              '.log',
              '.sh',
              '.txt',
              '.jpg',
              '.png',
              '.gif',
              '/fl_Input/',
              '/fl_Userinput',
              '/fl_Out_Std',
              '/fl_Plot',
              '/result.guess.lcao',
              '/result.guess.occ',
              '/guess.lcao',
              '/guess.occ',
              '/fl_Table/OrbitalTable',
              '/fl_Table/DensityTable',
              '/fl_Table/XcpotTable',
              '/fl_Work/fl_Mtr_Hpq.matrix',
              '/fl_Work/fl_Mtr_Spq.matrix',
              '/fl_Work/fl_Mtr_Sab.matrix',
              '/fl_Work/fl_Mtr_Sab2.matrix',
              '/fl_Work/fl_Mtr_Sgd.matrix',
              '/fl_Work/fl_Vct_Eigval',
              '/lo_Work/fl_Mtr_C.lo.',
              '/lo_Work/fl_Mtr_C.matrix.frag'
              )

alive_re_list = (re.compile('.*/fl_Work/fl_Mtr_C\.matrix\.\D+(\d+)'),
                 re.compile('.*/fl_Mtr_Fpq\.matrix\.\D+(\d+)')
                 )

def main():
    # parse options
    #
    parser = optparse.OptionParser(usage="%prog [directory]", version="%prog 1.0")
    # set option
    parser.add_option("-v", action="store_true", dest="verbose", default=False)

    (opts, args) = parser.parse_args()

    if (len(args) <= 0):
        parser.print_help()
        exit()
    elif (len(args) >= 1):
        pdf_dir = str(args[0])

    if (opts.verbose == True):
        print "top dir : %s" % (pdf_dir)

    #zip = zipfile.ZipFile('work-up.zip', 'w')
    #zip.create_system = 3 # for permission data
    arc = tarfile.TarFile.open('workup.tar.bz2', 'w:bz2')
    shrink_dir(pdf_dir, arc)
    for (root, dirs, files) in os.walk(pdf_dir):
        for d in dirs:
            path =  os.path.join(root, d)
            shrink_dir(path, arc)
    #zip.close()
    tarfile.close()


def shrink_dir(dir, arc):
    gbi = PdfGlobalinput(dir + "/fl_Input/fl_Globalinput")
    last_iter = gbi.get_value('SCF', 'control-iteration')
    if (last_iter != None):
        # do shrink!
        for search_path, dirs, files in os.walk(dir):
            for file in files:
                full_path = os.path.join(search_path, file)
                is_alive = check(full_path, last_iter)
                if (is_alive == False):
                    print "archived and deleted. %s" % (full_path)
                    #zip_object.write(full_path)
                    arc.add(full_path)
                    os.remove(full_path)

                
def check(full_path, last_iter):
    for item in alive_list:
        if (full_path.find(item) != -1):
            return True

    for re_item in alive_re_list:
        matchObj = re_item.match(full_path)
        if (matchObj != None):
            if (matchObj.group(1) == last_iter):
                return True

    return False


def make_tar(arcname, file, compression='bz2'):
    if compression:
        dest_ext = '.' + compression
    else:
        dest_ext = ''
    dest_name = '%s.tar%s' % (arcname, dest_ext)
    if compression:
        dest_cmp = ':' + compression
    else:
        dest_cmp = ''
    out = tarfile.TarFile.open(arcname, 'w'+dest_cmp)
    out.add(file)
    out.close()

if __name__ == '__main__':
    main()
    



