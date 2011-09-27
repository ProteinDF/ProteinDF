#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pack ProteinDF data
"""

import sqlite3


def SimpleMapper(object):
    """
    O/R mapper
    """

    
    @classmethod
    def get_connection(cls):
        return cls.connection
    


def main():
    con = sqlite3.connect(":memory")
    cur.execute("CREATE TABLE ")
    
if __name__ == '__main__':
    main()
