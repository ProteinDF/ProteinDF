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

#ifndef TLSYSTEM_H
#define TLSYSTEM_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <string>

class TlSystem {
 public:
  static int getPID();
  static int getPPID();
  static std::string getEnv(const std::string& key);

  /// ホスト名を返す
  static std::string getHostName();

  /// 使用された常駐セットサイズの最大値(MB)を返す
  static double getMaxRSS();

 private:
  static int pid_;
  static int ppid_;

  static const int MAX_HOSTNAME_LENGTH;
  static std::string hostname_;
};

#endif  // TLSYSTEM_H
