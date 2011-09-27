#include <cstdio>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>

#include "TlFile.h"

std::size_t TlFile::BUFFER_SIZE = 4096;

bool TlFile::isExist(const std::string& rFilePath)
{
    bool bAnswer = true;
    std::ifstream fs;
    fs.open(rFilePath.c_str());

    if (fs.is_open()) {
        bAnswer = true;
    } else {
        bAnswer = false;
    }

    fs.close();
    return bAnswer;
}


int TlFile::copy(const std::string& fromFilePath, const std::string& destFilePath)
{
    int answer = 0;

    if (TlFile::isExist(fromFilePath) == true) {
        std::ifstream in(fromFilePath.c_str(), std::ios_base::in | std::ios_base::binary);
        std::ofstream out(destFilePath.c_str(), std::ios_base::out | std::ios_base::binary);
        
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


int TlFile::remove(const std::string& filePath)
{
    int ans = 0;
    if (TlFile::isExist(filePath) == true) {
        ans = std::remove(filePath.c_str());
    }

    return ans;
}


int TlFile::rename(const std::string& oldName, const std::string& newName)
{
    int ans = 0;
    if (TlFile::isExist(oldName) == true) {
        ans = std::rename(oldName.c_str(), newName.c_str());
    }

    return ans;
}


size_t TlFile::getSize(const std::string& filePath)
{
    size_t answer = 0;

    struct stat fileInfo;
    if (stat(filePath.c_str(), &fileInfo) == 0) {
        answer = fileInfo.st_size;
    }

    return answer;
}
