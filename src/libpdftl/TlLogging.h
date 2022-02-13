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

#ifndef TLLOGGING_H
#define TLLOGGING_H

#include <fstream>
#include <string>

/// ログを扱うクラス
///
/// MPI並列プログラム用に出力レベルを2つ(master用, worker用)用意しています。
/// 出力レベル以上の出力が無ければ無視されます。
class TlLogging {
public:
    enum Level { TL_DEBUG,
                 TL_INFO,
                 TL_WARN,
                 TL_ERROR,
                 TL_CRITICAL };

public:
    static TlLogging& getInstance();

public:
    void setFilePath(const std::string& filePath);
    void setProcID(const int procID);

    TlLogging::Level getLevel(const bool isMaster = true) const {
        TlLogging::Level answer;
        if (isMaster) {
            answer = this->masterLevel_;
        } else {
            answer = this->workerLevel_;
        }
        return answer;
    }

    void setLevel(const TlLogging::Level masterLevel,
                  const TlLogging::Level workerLevel = TL_WARN);

    void debug(const std::string& msg);
    void info(const std::string& msg);
    void warn(const std::string& msg);
    void error(const std::string& msg);
    void critical(const std::string& msg);

private:
    TlLogging();
    TlLogging(const TlLogging& rhs);
    ~TlLogging();

    void setLevel();

    void open();
    void close();

    void output(const int level, const std::string& msg);
    std::string format(const int level, const std::string& input);

private:
    static TlLogging* pInstance_;
    static const char* pLevelStr_[];

    std::string filePath_;

    /// 出力用ファイルストリーム
    std::fstream output_;

    /// 出力レベル
    Level masterLevel_;
    Level workerLevel_;
    Level level_;

    /// プロセスIDを保持
    int procID_;
};

#endif  // TLLOGGING_H
