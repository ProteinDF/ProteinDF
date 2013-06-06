#ifndef PDFUSERINPUT_H
#define PDFUSERINPUT_H

#include <string>
#include "TlParameter.h"
#include "TlSerializeData.h"
#include "TlLogging.h"

/// 入力データの読み込み処理を行うクラス
class PdfUserInput {
public:
    PdfUserInput(const std::string& filePath = "fl_Userinput");
    ~PdfUserInput();

public:
    TlSerializeData getSerializeData() const;
    TlParameter getParameter() const;
    
    /// aliasの変換を行う
    void alias();
    bool check();

    void load();

public:
    TlSerializeData getBasisInfo(const std::string& basisName);
    
private:
    void load_conventional();

    void molecule_geometry_cartesian_input(const std::string& str);
    void moleculeBasisSetOrbital(const std::string& str);
    void moleculeBasisSetDensityAuxiliary(const std::string& str);
    void moleculeBasisSetExchangeAuxiliary(const std::string& str);
    void moleculeBasisSetGridFree(const std::string& str);

private:
    std::string filePath_;
    TlSerializeData data_;
    TlParameter param_;

    TlLogging& log_;
};

#endif // PDFUSERINPUT_H

