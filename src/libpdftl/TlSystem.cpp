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

#include "TlSystem.h"

#include <sys/mman.h>
#include <sys/types.h>

#include <cstdlib>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif  // HAVE_UNISTD_H

#ifdef HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif  // HAVE_SYS_RESOURCE_H

int TlSystem::pid_ = -1;
int TlSystem::ppid_ = -1;
const int TlSystem::MAX_HOSTNAME_LENGTH = 256;
std::string TlSystem::hostname_ = "";

int TlSystem::getPID() {
#ifdef HAVE_UNISTD_H
    if (TlSystem::pid_ == -1) {
        TlSystem::pid_ = getpid();
    }
#endif  // HAVE_UNISTD_H

    return TlSystem::pid_;
}

int TlSystem::getPPID() {
#ifdef HAVE_UNISTD_H
    if (TlSystem::ppid_ == -1) {
        TlSystem::ppid_ = getppid();
    }
#endif  // HAVE_UNISTD_H
    return TlSystem::ppid_;
}

std::string TlSystem::getEnv(const std::string& key) {
    std::string ans(std::getenv(key.c_str()));

    return ans;
}

std::string TlSystem::getHostName() {
#ifdef HAVE_UNISTD_H
    if (TlSystem::hostname_.empty()) {
        char* pHostName = new char[MAX_HOSTNAME_LENGTH];
        gethostname(pHostName, MAX_HOSTNAME_LENGTH);

        TlSystem::hostname_ = std::string(pHostName);

        delete[] pHostName;
        pHostName = NULL;
    }
#endif  // HAVE_UNISTD_H

    return TlSystem::hostname_;
}

double TlSystem::getMaxRSS() {
    double answer = 0.0;

#ifdef HAVE_SYS_RESOURCE_H
    struct rusage resourceUsage;
    getrusage(RUSAGE_SELF, &resourceUsage);
    answer = resourceUsage.ru_maxrss / 1024.0;
#endif  // HAVE_SYS_RESOURCE_H

    return answer;
}

char* TlSystem::newMmap(const std::size_t fileSize, const int fd) {
    // const std::size_t pageSize = sysconf(_SC_PAGE_SIZE);
    // const std::size_t mapSize = (this->fileSize_ / pageSize +1) * pageSize;

    const int prot = PROT_READ | PROT_WRITE;
    // const int flags = MAP_SHARED;  //  MAP_HUGETLB
    const int flags = MAP_PRIVATE;
    // const int flags = MAP_PRIVATE | MAP_HUGETLB;

    char* addr = reinterpret_cast<char*>(::mmap(NULL, fileSize, prot, flags, fd, 0));
    if (addr == MAP_FAILED) {
        do {
            perror("mmap");
            exit(EXIT_FAILURE);
        } while (0);
    }

    return addr;
}
