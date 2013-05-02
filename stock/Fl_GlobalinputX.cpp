#include "Fl_GlobalinputX.h"

////////////////////////////////////////////////////////////////////////
// Fl_GlobalinputX

Fl_GlobalinputX::Fl_GlobalinputX() : TlParameter()
{
    this->setFilePath("fl_Input/fl_Globalinput");
}

Fl_GlobalinputX::Fl_GlobalinputX(const Fl_GlobalinputX& rhs) : TlParameter(rhs)
{
}

Fl_GlobalinputX::~Fl_GlobalinputX()
{
}

bool Fl_GlobalinputX::load()
{
    return TlParameter::load("fl_Input/fl_Globalinput");
}

bool Fl_GlobalinputX::save() const
{
    return TlParameter::save("fl_Input/fl_Globalinput");
}


////////////////////////////////////////////////////////////////////////
// Fl_RestartX

Fl_RestartX::Fl_RestartX() : TlParameter()
{
    this->setFilePath("fl_Input/fl_Restart");

    (*this)["ProteinDF"]["DfIntegrals"       ] = "not done";
    (*this)["ProteinDF"]["DfInitialguess"    ] = "not done";
    (*this)["ProteinDF"]["DfScf"             ] = "not done";
    (*this)["ProteinDF"]["DfPhysics"         ] = "not done";
    (*this)["DfScf"]["dfScf-iteration-number"] = "1";
    (*this)["DfScf"]["write-norbcut"         ] = "not done";
}

Fl_RestartX::Fl_RestartX(const Fl_RestartX& rhs) : TlParameter(rhs)
{
}

Fl_RestartX::~Fl_RestartX()
{
}

bool Fl_RestartX::load()
{
    return TlParameter::load("fl_Input/fl_Restart");
}

bool Fl_RestartX::save() const
{
    return TlParameter::save("fl_Input/fl_Restart");
}

