#ifndef TLFILE_H
#define TLFILE_H

#include <string>

class TlFile {
public:
    /// rFilePathが存在していればtrue、なければfalseを返す。
    static bool isExist(const std::string& rFilePath);

    /// ファイルをコピーする。
    ///
    /// ファイルが存在しなかった場合は何もしない。
    static int copy(const std::string& fromFilePath, const std::string& destFilePath);

    /// 指定したパスのファイルを削除する。
    ///
    /// ファイルが存在しなかった場合は何もしない。
    /// @retval 0 成功
    /// @retval 0以外 エラー
    static int remove(const std::string& filePath);

    /// ファイル名を変更する。
    ///
    /// ファイルが存在しなかった場合は何もしない。
    /// @param oldName 変更前の名前
    /// @param newName 新しい名前
    /// @retval 0 成功
    /// @retval 0以外 エラー
    static int rename(const std::string& oldName, const std::string& newName);

    static size_t getSize(const std::string& filePath);

private:
    static std::size_t BUFFER_SIZE;
};

#endif // TLFILE_H
