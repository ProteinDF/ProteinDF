#ifndef FL_USERINPUTX_H
#define FL_USERINPUTX_H

#include <string>
#include "Fl_GlobalinputX.h"
#include "TlParameter.h"
#include "TlSerializeData.h"

/** 入力データの読み込み処理を行うクラス
 */
class Fl_UserinputX : public TlParameter {
public:
    Fl_UserinputX(const std::string& sFilePath = "fl_Userinput");
    ~Fl_UserinputX();

public:
    TlParameter getParameter() const;
    TlSerializeData getSerializeData() const;
    
    /// aliasの変換を行う
    void alias();
    bool check();
    void load();

    // void save();

private:

    std::string molecule_geometry_cartesian_input(const std::string& str);
    std::string moleculeBasisSetOrbital(const std::string& str);
    std::string moleculeBasisSetDensityAuxiliary(const std::string& str);
    std::string moleculeBasisSetExchangeAuxiliary(const std::string& str);

private:
    std::string m_sFilePath;
    TlParameter m_param;
    TlSerializeData data_;
};

#endif // FL_USERINPUTX_H

