// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <cassert>
#include <fstream>
#include <iostream>
#include <ostream>

// C
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "TlFile.h"

std::size_t TlFile::BUFFER_SIZE = 4096;

bool TlFile::isExistFile(const std::string& filePath) {
    bool answer = true;
    std::ifstream fs;
    fs.open(filePath.c_str());

    if (fs.is_open()) {
        answer = true;
    } else {
        answer = false;
    }

    fs.close();
    return answer;
}

int TlFile::copy(const std::string& fromFilePath,
                 const std::string& destFilePath) {
    int answer = 0;

    if (TlFile::isExistFile(fromFilePath) == true) {
        std::ifstream in(fromFilePath.c_str(),
                         std::ios_base::in | std::ios_base::binary);
        std::ofstream out(destFilePath.c_str(),
                          std::ios_base::out | std::ios_base::binary);

        char* buf = new char[TlFile::BUFFER_SIZE];

        do {
            in.read(buf, BUFFER_SIZE);
            out.write(buf, in.gcount());
        } while (in.gcount() > 0);

        in.close();
        out.close();

        delete[] buf;
        buf = NULL;
    }

    return answer;
}

int TlFile::remove(const std::string& filePath) {
    int ans = 0;
    if (TlFile::isExistFile(filePath) == true) {
        ans = std::remove(filePath.c_str());
    }

    return ans;
}

int TlFile::rename(const std::string& oldName, const std::string& newName) {
    int ans = 0;
    if (TlFile::isExistFile(newName) == true) {
        TlFile::remove(newName);
    }

    ans = std::rename(oldName.c_str(), newName.c_str());  // stdio.h
    if (ans != 0) {
        if (errno == EXDEV) {
            // both files exist different file system
            const int ret_cp = TlFile::copy(oldName, newName);
            if (ret_cp == 0) {
                int ret_rm = TlFile::remove(oldName);
                if (ret_rm == 0) {
                    ans = 0;
                }
            }
        }
    }

    return ans;
}

int TlFile::move(const std::string& oldName, const std::string& newName) {
    return TlFile::rename(oldName, newName);
}

std::size_t TlFile::getFileSize(const std::string& filePath) {
    size_t answer = 0;

    std::fstream fs;
    fs.open(filePath.c_str(), std::ios::in | std::ios::binary);
    fs.seekg(0, std::ios_base::end);
    std::streampos length = fs.tellg();
    if (length != -1) {
        answer = static_cast<std::size_t>(length);
    }
    fs.close();

    return answer;
}

std::string TlFile::getTempFilePath(const std::string& tmpdir) {
    const std::string baseName = tmpdir + std::string("pdftmp.XXXXXX");
    char* pFileName = new char[baseName.size() + 1];
    strncpy(pFileName, baseName.c_str(), baseName.size() + 1);

    int fd = mkstemp(pFileName);
    close(fd);

    const std::string fileName(pFileName);
    TlFile::remove(fileName);
    assert(TlFile::isExistFile(fileName) != true);

    delete[] pFileName;
    pFileName = NULL;

    return fileName;
}
