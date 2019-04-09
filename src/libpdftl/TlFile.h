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

#ifndef TLFILE_H
#define TLFILE_H

#include <cstdio>
#include <string>

class TlFile {
   public:
    /// filePathが存在していればtrue、なければfalseを返す。
    static bool isExistFile(const std::string& filePath);

    /// ファイルをコピーする。
    ///
    /// ファイルが存在しなかった場合は何もしない。
    static int copy(const std::string& fromFilePath,
                    const std::string& destFilePath);

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

    /// get file size (bytes).
    static std::size_t getFileSize(const std::string& filePath);

    /// return temporary file path
    static std::string getTempFilePath(const std::string& tmpdir = "/tmp/");

   private:
    static std::size_t BUFFER_SIZE;
};

#endif  // TLFILE_H
