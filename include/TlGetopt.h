#ifndef TLGETOPT_H
#define TLGETOPT_H

#include <map>
#include <string>

/**
 *  コマンドラインオプションを処理するクラス
 *
 *  コンストラクタにargc, argvを設定し、listに所定の書式で
 *  読み込むオプションを指定すると、このオブジェクトに対し、
 *  []メソッドでオプションの値にアクセスできる。
 *  オプションの書式は、いわゆるgetoptの書式に準ずる。
 *  long optionは指定できない。
 *  また、引数の値には、0から始まるインデックスで、
 *  同じく[]メソッドでアクセスできる。
 *  argv[0]の値はコマンド名になる前提で、
 *  オプション解析をargv[1]から処理しているが、
 *  argv[0]がコマンド名になるかは、実装定義であるので注意が必要である。
 */
class TlGetopt {
public:
    TlGetopt(int argc, char* argv[], const char* list);
    ~TlGetopt();

    // accession
public:
    int getCount() const {
        return this->m_nCount;
    };
    const std::string operator[](const std::string& sKey) const;
    const std::string operator[](unsigned int n) const;

protected:
    void initialize();
    void parseArgv(int argc, char* argv[], const char* list);

protected:
    /**
     *  引数の数
     */
    int m_nCount;

    /**
     *  データ保持用
     */
    std::map<std::string, std::string> m_Data;

    /**
     *  エラー文字列
     */
    std::string m_sError;
};

#endif // TLGETOPT_H
