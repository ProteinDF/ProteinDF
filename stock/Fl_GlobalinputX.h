#ifndef FL_GLOBALINPUTX_H
#define FL_GLOBALINPUTX_H

#include "TlParameter.h"

////////////////////////////////////////////////////////////////////////
// class Fl_GlobalinputX
//

class Fl_GlobalinputX : public TlParameter {
public:
    Fl_GlobalinputX();
    Fl_GlobalinputX(const Fl_GlobalinputX& rhs);
    virtual ~Fl_GlobalinputX();

public:
    bool load();
    bool save() const;
};


////////////////////////////////////////////////////////////////////////
// class Fl_RestartX
//
class Fl_RestartX : public TlParameter {
public:
    explicit Fl_RestartX();
    Fl_RestartX(const Fl_RestartX& rhs);
    virtual ~Fl_RestartX();

public:
    bool load();
    bool save() const;

private:
    bool m_bAutoLoad;
    bool m_bAutoSave;
};

#endif // FL_GLOBALINPUTX_H
