#!/usr/bin/env python

import getopt, sys
import TlMatrix

def main():
    try:
        optlist, args = getopt.gnu_getopt(sys.argv[1:], "hs:", longopts=["help", "size="])
    except getopt.GetoptError:
        # exit before display the help message
        usage()
        sys.exit(2)

    for opt, arg in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        if opt in ("-s", "--size"):
            print "size = ", arg

    if (len(args) == 0):
        usage()
        sys.exit(2)
        
    file_path1 = args.pop(0)
    print "reading matrix file: ", file_path1

    matrix_box1 = TlMatrix.TlMatrix()
    matrix_sym1 = TlMatrix.TlMatrix_Symmetric()

    if (matrix_box1.is_loadable(file_path1) == True):
        matrix_box1.load_binary(file_path1)
    elif (matrix_sym1.is_loadable(file_path1) == True):
        matrix_sym1.load_binary(file_path1)
    else:
        print "unkown file type. stop."
    
    file_path2 = args.pop(0)
    print "reading matrix file: ", file_path2

    matrix_box2 = TlMatrix.TlMatrix()
    matrix_sym2 = TlMatrix.TlMatrix_Symmetric()

    if (matrix_box2.is_loadable(file_path2) == True):
        matrix_box2.load_binary(file_path2)
    elif (matrix_sym2.is_loadable(file_path2) == True):
        matrix_sym2.load_binary(file_path2)
    else:
        print "unkown file type. stop."

    matrix_sym3 = TlMatrix.TlMatrix_Symmetric()
    matrix_sym3 = matrix_sym1 - matrix_sym2
    matrix_sym3.printout()


def usage():
    print "to implement!"


if __name__ == '__main__':
    main()
    
