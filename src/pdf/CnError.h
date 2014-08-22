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

#ifndef CNERROR_H
#define CNERROR_H

#include <string>

/// エラー時のメッセージと終了処理を扱うクラス
class CnError {
public:
    /// デフォルトコンストラクタ
    CnError();

    /// デストラクタ
    ~CnError();

public:
    /// 処理を中止させる
    void abort(const std::string& sMsg = "");

    // to terminate program, abort of Standard C Library is called in following member function
    void abort(const std::string& ClassName, const std::string& ObjName, const std::string& MemFunc, const std::string& Message);
};

extern CnError CnErr;   // global object of CnError Class

#endif // CNERROR_H

