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

#ifndef FL_FRAGMENT_H
#define FL_FRAGMENT_H

#include <map>
#include <string>
#include <vector>
#include "Fl_Geometry.h"

class Fl_Fragment {
   public:
    Fl_Fragment(const Fl_Geometry& flGeom);
    ~Fl_Fragment();

   public:
    int getNumOfFragments() const;
    int getNumOfOrbitals(int fragindex) const;
    int getNumOfOccupiedOrbitals(int fragindex) const;
    int getNumOfOccupiedAlphaOrbitals(int fragindex) const;
    int getNumOfOccupiedBetaOrbitals(int fragindex) const;

    void calcDefault(const std::vector<int>& norbcut);
    // int calcDefault(int* norbcut);

   private:
    std::string getStr() const;
    void save() const;
    void load();

   private:
    const Fl_Geometry& flGeom_;
    //   void  level(int pout);
    //   int   read();
    //   void  show() const;
    //   void  show(const std::string& title) const;

    // put or set data of fragment.
    //   int       putFramename(std::string& name);
    //   int       putNumberFragment(int value);
    //   int       putCountFragment(int value);
    //   int       putIndex(int fragnum, int value);
    //   int       putName(int fragnum, std::string& name);
    //   int       putNumberElectron(int fragnum, int value);
    //   int       putNumberAlphaelectron(int fragnum, int value);
    //   int       putNumberBetaelectron(int fragnum, int value);
    //   int       putNumberOrbital(int fragnum, int value);
    //   int       putNumberOccupied(int fragnum, int value);
    //   int       putNumberOccupiedAlpha(int fragnum, int value);
    //   int       putNumberOccupiedBeta(int fragnum, int value);
    //   int       putNumberUnoccupied(int fragnum, int value);
    //   int       putNumberUnoccupiedAlpha(int fragnum, int value);
    //   int       putNumberUnoccupiedBeta(int fragnum, int value);

    // get data of fragment.
    //   std::string& getFramename();
    //   int      getCountFragment();
    //   std::string& getName(int fragindex);
    //   int      getNumberElectron(int fragindex);
    //   int      getNumberAlphaelectron(int fragindex);
    //   int      getNumberBetaelectron(int fragindex);
    //   int      getNumberUnoccupied(int fragindex);
    //   int      getNumberUnoccupiedAlpha(int fragindex);
    //   int      getNumberUnoccupiedBeta(int fragindex);

    //   void  changeNumberFragment(int value);

   private:
    //   enum { ON=1, OFF=0, DEBUG=OFF };
    //   enum { ExitValue = 1 };  // Parameter of Standard C Library exit().

    std::string m_filePath;

    struct Fragment {
       public:
        Fragment()
            : numOfElectrons(-1),
              numOfAlphaElectrons(-1),
              numOfBetaElectrons(-1),
              numOfOrbitals(-1),
              numOfOccupiedOrbitals(-1),
              numOfOccupiedAlphaOrbitals(-1),
              numOfOccupiedBetaOrbitals(-1),
              numOfUnoccupiedOrbitals(-1),
              numOfUnoccupiedAlphaOrbitals(-1),
              numOfUnoccupiedBetaOrbitals(-1) {}

       public:
        // int index;                   // index of the fragment
        std::string name;                // name of the fragment
        int numOfElectrons;              // number of electrons       (for RKS)
        int numOfAlphaElectrons;         // number of alpha electrons (for UKS)
        int numOfBetaElectrons;          // number of beta electrons  (for UKS)
        int numOfOrbitals;               // number of orbitals        (for RKS)
        int numOfOccupiedOrbitals;       // number of occupied orbitals
        int numOfOccupiedAlphaOrbitals;  // number of occupied orbitals
        int numOfOccupiedBetaOrbitals;   // number of occupied orbitals
        int numOfUnoccupiedOrbitals;     // number of unoccupied orbitals
        int numOfUnoccupiedAlphaOrbitals;  // number of unoccupied orbitals
        int numOfUnoccupiedBetaOrbitals;   // number of unoccupied orbitals
    };

    std::map<int, Fragment> m_fragments;
    std::string m_frameName;  // name of frame molecule

    //   int*     index_fragnum;
    //   int      number_fragment;       // maximum number of fragments
    //   int      count_fragment;        // number of registered fragments
};

#endif  // FL_FLAGMENT_H
