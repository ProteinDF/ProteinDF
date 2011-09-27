#ifndef FL_VCT_OTABLEX_H
#define FL_VCT_OTABLEX_H

#include <vector>
#include <string>

class Fl_Vct_OtableX {
public:
    struct Info {
public:
        explicit Info(int nAtomBelong =0, int nBasisType =0, int nContract =0, double dPreExp =0.0, double dExpAlpha =0.0) :
                atomBelong(nAtomBelong), basisType(nBasisType), contract(nContract), preExp(dPreExp), expAlpha(dExpAlpha) {
        }

        int atomBelong;
        int basisType;
        int contract;
        double preExp;
        double expAlpha;
    };

public:
    Fl_Vct_OtableX();
    ~Fl_Vct_OtableX();

public:
    int getSize() const;

    void resize(int size);
    void addData(const Info& info);

    Info  operator[](int index) const;
    Info& operator[](int index);

public:
    bool save() const;
    bool load();

protected:
    bool save(std::ofstream& ofs) const;
    bool load(std::ifstream& ifs);

private:
    std::string m_sFilePath;
    std::vector<Info> m_data;
};

#endif // FL_VCT_OTABLEX_H
