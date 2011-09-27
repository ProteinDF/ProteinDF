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

