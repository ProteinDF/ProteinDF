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

#include "DfGrid.h"

#define SQ2 1.414213562373095049
#define SQ1_2 0.707106781186547524
#define SQ1_3 0.577350269189625765
#define SMALL 1
#define TOOBIG 30.0
#define EPS 0.000000001 /* original #define EPS         0.000000000001 */

void points(int grid, int nOgrid, double r0, double* Ogridx, double* Ogridy,
            double* Ogridz, double* w, double weight) {
    // Set the grid coordinates in "Ogridr" and the Lebedev weight in "w"
    // Get atom core coordinates

    double corex, corey, corez;
    int i, j;

    double A1, A2, A3;
    double B1, B2, B3, B4;
    double C1;
    double D1;

    double p1;
    double q1;

    double bm[4];
    double bl[4];
    double dr;
    double du;
    double dw;

    double CORRP = 1.00;

    corex = Ogridx[0];
    corey = Ogridy[0];
    corez = Ogridz[0];

    for (i = 0; i < nOgrid; ++i) {
        w[i] = weight;
        Ogridx[i] = Ogridy[i] = Ogridz[i] = 0.0;
    }

    if (grid == -1 || grid == 2) {
        cout << "Sorry, do not support now." << endl;

    } else if (grid >= 3) {
        //  SG-1 Metod
        //*******Set the Lebedev Grid of 6 points**********************

        if (nOgrid == 6) {
            // Set the Lebedev weight
            for (i = 0; i < nOgrid; ++i) w[i] *= 0.166666666666667;

            // Set the grid coordinates
            Ogridx[0] = 1.0;
            Ogridy[0] = 0.0;
            Ogridz[0] = 0.0;
            Ogridx[1] = 0.0;
            Ogridy[1] = 1.0;
            Ogridz[1] = 0.0;
            Ogridx[2] = 0.0;
            Ogridy[2] = 0.0;
            Ogridz[2] = 1.0;
            Ogridx[3] = -1.0;
            Ogridy[3] = 0.0;
            Ogridz[3] = 0.0;
            Ogridx[4] = 0.0;
            Ogridy[4] = -1.0;
            Ogridz[4] = 0.0;
            Ogridx[5] = 0.0;
            Ogridy[5] = 0.0;
            Ogridz[5] = -1.0;

            // FOR CHECK
            //      for ( i = 0; i < nOgrid; ++i ) {
            //        cout << i  << " " <<  Ogridx[i] << "  " << Ogridy[i] << "
            //        " << Ogridz[i] << "  "
            //             << (
            //             Ogridx[i]*Ogridx[i]+Ogridy[i]*Ogridy[i]+Ogridz[i]*Ogridz[i]
            //             ) << "  " << w[i] <<  endl;
            //      }
            // END
            // Scale and Shift
            for (i = 0; i < nOgrid; ++i) {
                Ogridx[i] = r0 * Ogridx[i] + corex;
                Ogridy[i] = r0 * Ogridy[i] + corey;
                Ogridz[i] = r0 * Ogridz[i] + corez;
            }

            //******Set the Lebedev Grid of 38 points*******************

        } else if (nOgrid == 38) {
            // Set the Lebedev parameter

            A1 = 0.00952380952387;
            A3 = 0.0321428571429;
            C1 = 0.0285714285714;
            p1 = 0.888073833977;
            q1 = 0.459700843381;

            // Set the Lebedev weight
            for (i = 0; i < 6; ++i) w[i] *= A1;

            for (i = 6; i < 14; ++i) w[i] *= A3;

            for (i = 14; i < 38; ++i) w[i] *= C1;

            // Set the grid coordinates
            Ogridx[0] = 1.0;
            Ogridy[0] = 0.0;
            Ogridz[0] = 0.0;
            Ogridx[1] = 0.0;
            Ogridy[1] = 1.0;
            Ogridz[1] = 0.0;
            Ogridx[2] = 0.0;
            Ogridy[2] = 0.0;
            Ogridz[2] = 1.0;
            Ogridx[3] = -1.0;
            Ogridy[3] = 0.0;
            Ogridz[3] = 0.0;
            Ogridx[4] = 0.0;
            Ogridy[4] = -1.0;
            Ogridz[4] = 0.0;
            Ogridx[5] = 0.0;
            Ogridy[5] = 0.0;
            Ogridz[5] = -1.0;

            Ogridx[6] = SQ1_3;
            Ogridy[6] = SQ1_3;
            Ogridz[6] = SQ1_3;
            Ogridx[7] = -SQ1_3;
            Ogridy[7] = SQ1_3;
            Ogridz[7] = SQ1_3;
            Ogridx[8] = SQ1_3;
            Ogridy[8] = -SQ1_3;
            Ogridz[8] = SQ1_3;
            Ogridx[9] = SQ1_3;
            Ogridy[9] = SQ1_3;
            Ogridz[9] = -SQ1_3;
            Ogridx[10] = -SQ1_3;
            Ogridy[10] = -SQ1_3;
            Ogridz[10] = SQ1_3;
            Ogridx[11] = -SQ1_3;
            Ogridy[11] = SQ1_3;
            Ogridz[11] = -SQ1_3;
            Ogridx[12] = SQ1_3;
            Ogridy[12] = -SQ1_3;
            Ogridz[12] = -SQ1_3;
            Ogridx[13] = -SQ1_3;
            Ogridy[13] = -SQ1_3;
            Ogridz[13] = -SQ1_3;

            // set ci1
            Ogridx[14] = p1;
            Ogridy[14] = q1;
            Ogridz[14] = 0.0;
            Ogridx[15] = p1;
            Ogridy[15] = -q1;
            Ogridz[15] = 0.0;
            Ogridx[16] = -p1;
            Ogridy[16] = q1;
            Ogridz[16] = 0.0;
            Ogridx[17] = -p1;
            Ogridy[17] = -q1;
            Ogridz[17] = 0.0;

            Ogridx[18] = p1;
            Ogridy[18] = 0.0;
            Ogridz[18] = q1;
            Ogridx[19] = p1;
            Ogridy[19] = 0.0;
            Ogridz[19] = -q1;
            Ogridx[20] = -p1;
            Ogridy[20] = 0.0;
            Ogridz[20] = q1;
            Ogridx[21] = -p1;
            Ogridy[21] = 0.0;
            Ogridz[21] = -q1;

            Ogridx[22] = 0.0;
            Ogridy[22] = p1;
            Ogridz[22] = q1;
            Ogridx[23] = 0.0;
            Ogridy[23] = p1;
            Ogridz[23] = -q1;
            Ogridx[24] = 0.0;
            Ogridy[24] = -p1;
            Ogridz[24] = q1;
            Ogridx[25] = 0.0;
            Ogridy[25] = -p1;
            Ogridz[25] = -q1;

            Ogridx[26] = q1;
            Ogridy[26] = p1;
            Ogridz[26] = 0.0;
            Ogridx[27] = q1;
            Ogridy[27] = -p1;
            Ogridz[27] = 0.0;
            Ogridx[28] = -q1;
            Ogridy[28] = p1;
            Ogridz[28] = 0.0;
            Ogridx[29] = -q1;
            Ogridy[29] = -p1;
            Ogridz[29] = 0.0;

            Ogridx[30] = q1;
            Ogridy[30] = 0.0;
            Ogridz[30] = p1;
            Ogridx[31] = q1;
            Ogridy[31] = 0.0;
            Ogridz[31] = -p1;
            Ogridx[32] = -q1;
            Ogridy[32] = 0.0;
            Ogridz[32] = p1;
            Ogridx[33] = -q1;
            Ogridy[33] = 0.0;
            Ogridz[33] = -p1;

            Ogridx[34] = 0.0;
            Ogridy[34] = q1;
            Ogridz[34] = p1;
            Ogridx[35] = 0.0;
            Ogridy[35] = q1;
            Ogridz[35] = -p1;
            Ogridx[36] = 0.0;
            Ogridy[36] = -q1;
            Ogridz[36] = p1;
            Ogridx[37] = 0.0;
            Ogridy[37] = -q1;
            Ogridz[37] = -p1;

            // FOR CHECK
            //      for ( i = 0; i < nOgrid; ++i ) {
            //          cout << i << " " <<  Ogridx[i] << "  " << Ogridy[i] << "
            //          " << Ogridz[i] << "  "
            //               << (
            //               Ogridx[i]*Ogridx[i]+Ogridy[i]*Ogridy[i]+Ogridz[i]*Ogridz[i]
            //               ) << "  " << w[i] << endl;
            //      }
            // END
            // Scale and Shift
            for (i = 0; i < nOgrid; ++i) {
                Ogridx[i] = r0 * Ogridx[i] + corex;
                Ogridy[i] = r0 * Ogridy[i] + corey;
                Ogridz[i] = r0 * Ogridz[i] + corez;
            }

            // ******Set the Lebedev grid of 86 points*****************
        } else if (nOgrid == 86) {
            // Set the Lebedev parameter

            A1 = 0.0115440115441;
            A3 = 0.0119439090859;
            B1 = 0.0111105557106;
            B2 = 0.0118765012945;
            C1 = 0.0118123037469;

            bm[0] = 0.852518311701;
            bm[1] = 0.189063552885;

            bl[0] = 0.369602846454;
            bl[1] = 0.694354006603;

            p1 = 0.927330657151;
            q1 = 0.374243039090;

            // Set the Lebedev weight
            for (i = 0; i < 6; ++i) w[i] *= A1;

            for (i = 6; i < 14; ++i) w[i] *= A3;

            for (i = 14; i < 38; ++i) w[i] *= B1;

            for (i = 38; i < 62; ++i) w[i] *= B2;

            for (i = 62; i < 86; ++i) w[i] *= C1;

            // Set the grid coordinates
            Ogridx[0] = 1.0;
            Ogridy[0] = 0.0;
            Ogridz[0] = 0.0;
            Ogridx[1] = 0.0;
            Ogridy[1] = 1.0;
            Ogridz[1] = 0.0;
            Ogridx[2] = 0.0;
            Ogridy[2] = 0.0;
            Ogridz[2] = 1.0;
            Ogridx[3] = -1.0;
            Ogridy[3] = 0.0;
            Ogridz[3] = 0.0;
            Ogridx[4] = 0.0;
            Ogridy[4] = -1.0;
            Ogridz[4] = 0.0;
            Ogridx[5] = 0.0;
            Ogridy[5] = 0.0;
            Ogridz[5] = -1.0;

            Ogridx[6] = SQ1_3;
            Ogridy[6] = SQ1_3;
            Ogridz[6] = SQ1_3;
            Ogridx[7] = -SQ1_3;
            Ogridy[7] = SQ1_3;
            Ogridz[7] = SQ1_3;
            Ogridx[8] = SQ1_3;
            Ogridy[8] = -SQ1_3;
            Ogridz[8] = SQ1_3;
            Ogridx[9] = SQ1_3;
            Ogridy[9] = SQ1_3;
            Ogridz[9] = -SQ1_3;
            Ogridx[10] = -SQ1_3;
            Ogridy[10] = -SQ1_3;
            Ogridz[10] = SQ1_3;
            Ogridx[11] = -SQ1_3;
            Ogridy[11] = SQ1_3;
            Ogridz[11] = -SQ1_3;
            Ogridx[12] = SQ1_3;
            Ogridy[12] = -SQ1_3;
            Ogridz[12] = -SQ1_3;
            Ogridx[13] = -SQ1_3;
            Ogridy[13] = -SQ1_3;
            Ogridz[13] = -SQ1_3;

            // set bik
            for (i = 0; i < 2; ++i) {
                // bl[i] = sqrt( 1.0 - bm[i] * bm[i] ) / SQ2;

                Ogridx[14 + 24 * i + 0] = bl[i];
                Ogridy[14 + 24 * i + 0] = bl[i];
                Ogridz[14 + 24 * i + 0] = bm[i];
                Ogridx[14 + 24 * i + 1] = -bl[i];
                Ogridy[14 + 24 * i + 1] = bl[i];
                Ogridz[14 + 24 * i + 1] = bm[i];
                Ogridx[14 + 24 * i + 2] = bl[i];
                Ogridy[14 + 24 * i + 2] = -bl[i];
                Ogridz[14 + 24 * i + 2] = bm[i];
                Ogridx[14 + 24 * i + 3] = bl[i];
                Ogridy[14 + 24 * i + 3] = bl[i];
                Ogridz[14 + 24 * i + 3] = -bm[i];
                Ogridx[14 + 24 * i + 4] = -bl[i];
                Ogridy[14 + 24 * i + 4] = -bl[i];
                Ogridz[14 + 24 * i + 4] = bm[i];
                Ogridx[14 + 24 * i + 5] = -bl[i];
                Ogridy[14 + 24 * i + 5] = bl[i];
                Ogridz[14 + 24 * i + 5] = -bm[i];
                Ogridx[14 + 24 * i + 6] = bl[i];
                Ogridy[14 + 24 * i + 6] = -bl[i];
                Ogridz[14 + 24 * i + 6] = -bm[i];
                Ogridx[14 + 24 * i + 7] = -bl[i];
                Ogridy[14 + 24 * i + 7] = -bl[i];
                Ogridz[14 + 24 * i + 7] = -bm[i];
            }
            for (i = 0; i < 2; ++i) {
                for (j = 0; j < 8; ++j) {
                    Ogridx[14 + 24 * i + 8 + j] = Ogridy[14 + 24 * i + j];
                    Ogridy[14 + 24 * i + 8 + j] = Ogridz[14 + 24 * i + j];
                    Ogridz[14 + 24 * i + 8 + j] = Ogridx[14 + 24 * i + j];

                    Ogridx[14 + 24 * i + 16 + j] = Ogridz[14 + 24 * i + j];
                    Ogridy[14 + 24 * i + 16 + j] = Ogridx[14 + 24 * i + j];
                    Ogridz[14 + 24 * i + 16 + j] = Ogridy[14 + 24 * i + j];
                }
            }

            // set ci1
            Ogridx[62] = p1;
            Ogridy[62] = q1;
            Ogridz[62] = 0.0;
            Ogridx[63] = p1;
            Ogridy[63] = -q1;
            Ogridz[63] = 0.0;
            Ogridx[64] = -p1;
            Ogridy[64] = q1;
            Ogridz[64] = 0.0;
            Ogridx[65] = -p1;
            Ogridy[65] = -q1;
            Ogridz[65] = 0.0;

            Ogridx[66] = p1;
            Ogridy[66] = 0.0;
            Ogridz[66] = q1;
            Ogridx[67] = p1;
            Ogridy[67] = 0.0;
            Ogridz[67] = -q1;
            Ogridx[68] = -p1;
            Ogridy[68] = 0.0;
            Ogridz[68] = q1;
            Ogridx[69] = -p1;
            Ogridy[69] = 0.0;
            Ogridz[69] = -q1;

            Ogridx[70] = 0.0;
            Ogridy[70] = p1;
            Ogridz[70] = q1;
            Ogridx[71] = 0.0;
            Ogridy[71] = p1;
            Ogridz[71] = -q1;
            Ogridx[72] = 0.0;
            Ogridy[72] = -p1;
            Ogridz[72] = q1;
            Ogridx[73] = 0.0;
            Ogridy[73] = -p1;
            Ogridz[73] = -q1;

            Ogridx[74] = q1;
            Ogridy[74] = p1;
            Ogridz[74] = 0.0;
            Ogridx[75] = q1;
            Ogridy[75] = -p1;
            Ogridz[75] = 0.0;
            Ogridx[76] = -q1;
            Ogridy[76] = p1;
            Ogridz[76] = 0.0;
            Ogridx[77] = -q1;
            Ogridy[77] = -p1;
            Ogridz[77] = 0.0;

            Ogridx[78] = q1;
            Ogridy[78] = 0.0;
            Ogridz[78] = p1;
            Ogridx[79] = q1;
            Ogridy[79] = 0.0;
            Ogridz[79] = -p1;
            Ogridx[80] = -q1;
            Ogridy[80] = 0.0;
            Ogridz[80] = p1;
            Ogridx[81] = -q1;
            Ogridy[81] = 0.0;
            Ogridz[81] = -p1;

            Ogridx[82] = 0.0;
            Ogridy[82] = q1;
            Ogridz[82] = p1;
            Ogridx[83] = 0.0;
            Ogridy[83] = q1;
            Ogridz[83] = -p1;
            Ogridx[84] = 0.0;
            Ogridy[84] = -q1;
            Ogridz[84] = p1;
            Ogridx[85] = 0.0;
            Ogridy[85] = -q1;
            Ogridz[85] = -p1;

            // FOR CHECK
            //      for ( i = 0; i < nOgrid; ++i ) {
            //              cout <<  i  << " " <<  Ogridx[i] << "  " <<
            //              Ogridy[i] << " " << Ogridz[i] << "  "
            //                   << (
            //                   Ogridx[i]*Ogridx[i]+Ogridy[i]*Ogridy[i]+Ogridz[i]*Ogridz[i]
            //                   ) << "  " << w[i] << endl;
            //      }
            // END
            // Scale and Shift
            for (i = 0; i < nOgrid; ++i) {
                Ogridx[i] = r0 * Ogridx[i] + corex;
                Ogridy[i] = r0 * Ogridy[i] + corey;
                Ogridz[i] = r0 * Ogridz[i] + corez;
            }

            // ******Set the Lebedev grid of 194 points******************
        } else if (nOgrid == 194) {
            // Set the Lebedev parameter
            A1 = 0.00178234044724;
            A2 = 0.00571690594998;
            A3 = 0.00557338317884;
            B1 = 0.00551877146727;
            B2 = 0.00515823771181;
            B3 = 0.00560870408259;
            B4 = 0.00410677702817;
            C1 = 0.00505184606462;
            D1 = 0.00553024891623;

            bm[0] = 0.777493219315;
            bm[1] = 0.912509096867;
            bm[2] = 0.314196994183;
            bm[3] = 0.982972302707;

            bl[0] = 0.444693317871;
            bl[1] = 0.289246562758;
            bl[2] = 0.671297344270;
            bl[3] = 0.129933544765;

            p1 = 0.938319218138;
            q1 = 0.345770219761;

            dr = 0.836036015482;
            du = 0.159041710538;
            dw = 0.525118572443;

            // Set the Lebedev weight

            for (i = 0; i < 6; ++i) w[i] *= A1;

            for (i = 6; i < 18; ++i) w[i] *= A2;

            for (i = 18; i < 26; ++i) w[i] *= A3;

            for (i = 26; i < 50; ++i) w[i] *= B1;

            for (i = 50; i < 74; ++i) w[i] *= B2;

            for (i = 74; i < 98; ++i) w[i] *= B3;

            for (i = 98; i < 122; ++i) w[i] *= B4;

            for (i = 122; i < 146; ++i) w[i] *= C1;

            for (i = 146; i < 194; ++i) w[i] *= D1;

            // Set the grid coordinates
            Ogridx[0] = 1.0;
            Ogridy[0] = 0.0;
            Ogridz[0] = 0.0;
            Ogridx[1] = 0.0;
            Ogridy[1] = 1.0;
            Ogridz[1] = 0.0;
            Ogridx[2] = 0.0;
            Ogridy[2] = 0.0;
            Ogridz[2] = 1.0;
            Ogridx[3] = -1.0;
            Ogridy[3] = 0.0;
            Ogridz[3] = 0.0;
            Ogridx[4] = 0.0;
            Ogridy[4] = -1.0;
            Ogridz[4] = 0.0;
            Ogridx[5] = 0.0;
            Ogridy[5] = 0.0;
            Ogridz[5] = -1.0;

            Ogridx[6] = SQ1_2;
            Ogridy[6] = SQ1_2;
            Ogridz[6] = 0.0;
            Ogridx[7] = -SQ1_2;
            Ogridy[7] = SQ1_2;
            Ogridz[7] = 0.0;
            Ogridx[8] = SQ1_2;
            Ogridy[8] = -SQ1_2;
            Ogridz[8] = 0.0;
            Ogridx[9] = -SQ1_2;
            Ogridy[9] = -SQ1_2;
            Ogridz[9] = 0.0;

            for (i = 0; i < 4; ++i) {
                Ogridx[6 + 4 + i] = Ogridy[6 + i];
                Ogridy[6 + 4 + i] = Ogridz[6 + i];
                Ogridz[6 + 4 + i] = Ogridx[6 + i];

                Ogridx[6 + 8 + i] = Ogridz[6 + i];
                Ogridy[6 + 8 + i] = Ogridx[6 + i];
                Ogridz[6 + 8 + i] = Ogridy[6 + i];
            }

            Ogridx[18] = SQ1_3;
            Ogridy[18] = SQ1_3;
            Ogridz[18] = SQ1_3;
            Ogridx[19] = -SQ1_3;
            Ogridy[19] = SQ1_3;
            Ogridz[19] = SQ1_3;
            Ogridx[20] = SQ1_3;
            Ogridy[20] = -SQ1_3;
            Ogridz[20] = SQ1_3;
            Ogridx[21] = SQ1_3;
            Ogridy[21] = SQ1_3;
            Ogridz[21] = -SQ1_3;
            Ogridx[22] = -SQ1_3;
            Ogridy[22] = -SQ1_3;
            Ogridz[22] = SQ1_3;
            Ogridx[23] = -SQ1_3;
            Ogridy[23] = SQ1_3;
            Ogridz[23] = -SQ1_3;
            Ogridx[24] = SQ1_3;
            Ogridy[24] = -SQ1_3;
            Ogridz[24] = -SQ1_3;
            Ogridx[25] = -SQ1_3;
            Ogridy[25] = -SQ1_3;
            Ogridz[25] = -SQ1_3;

            // set bik; no need to set cik
            for (i = 0; i < 4; ++i) {
                // bl[i] = sqrt( 1.0 - bm[i] * bm[i] ) / SQ2;

                Ogridx[26 + 24 * i + 0] = bl[i];
                Ogridy[26 + 24 * i + 0] = bl[i];
                Ogridz[26 + 24 * i + 0] = bm[i];
                Ogridx[26 + 24 * i + 1] = -bl[i];
                Ogridy[26 + 24 * i + 1] = bl[i];
                Ogridz[26 + 24 * i + 1] = bm[i];
                Ogridx[26 + 24 * i + 2] = bl[i];
                Ogridy[26 + 24 * i + 2] = -bl[i];
                Ogridz[26 + 24 * i + 2] = bm[i];
                Ogridx[26 + 24 * i + 3] = bl[i];
                Ogridy[26 + 24 * i + 3] = bl[i];
                Ogridz[26 + 24 * i + 3] = -bm[i];
                Ogridx[26 + 24 * i + 4] = -bl[i];
                Ogridy[26 + 24 * i + 4] = -bl[i];
                Ogridz[26 + 24 * i + 4] = bm[i];
                Ogridx[26 + 24 * i + 5] = -bl[i];
                Ogridy[26 + 24 * i + 5] = bl[i];
                Ogridz[26 + 24 * i + 5] = -bm[i];
                Ogridx[26 + 24 * i + 6] = bl[i];
                Ogridy[26 + 24 * i + 6] = -bl[i];
                Ogridz[26 + 24 * i + 6] = -bm[i];
                Ogridx[26 + 24 * i + 7] = -bl[i];
                Ogridy[26 + 24 * i + 7] = -bl[i];
                Ogridz[26 + 24 * i + 7] = -bm[i];
            }

            for (i = 0; i < 4; ++i) {
                for (j = 0; j < 8; ++j) {
                    Ogridx[26 + 24 * i + 8 + j] = Ogridy[26 + 24 * i + j];
                    Ogridy[26 + 24 * i + 8 + j] = Ogridz[26 + 24 * i + j];
                    Ogridz[26 + 24 * i + 8 + j] = Ogridx[26 + 24 * i + j];

                    Ogridx[26 + 24 * i + 16 + j] = Ogridz[26 + 24 * i + j];
                    Ogridy[26 + 24 * i + 16 + j] = Ogridx[26 + 24 * i + j];
                    Ogridz[26 + 24 * i + 16 + j] = Ogridy[26 + 24 * i + j];
                }
            }

            // set ci1
            Ogridx[122] = p1;
            Ogridy[122] = q1;
            Ogridz[122] = 0.0;
            Ogridx[123] = p1;
            Ogridy[123] = -q1;
            Ogridz[123] = 0.0;
            Ogridx[124] = -p1;
            Ogridy[124] = q1;
            Ogridz[124] = 0.0;
            Ogridx[125] = -p1;
            Ogridy[125] = -q1;
            Ogridz[125] = 0.0;

            Ogridx[126] = p1;
            Ogridy[126] = 0.0;
            Ogridz[126] = q1;
            Ogridx[127] = p1;
            Ogridy[127] = 0.0;
            Ogridz[127] = -q1;
            Ogridx[128] = -p1;
            Ogridy[128] = 0.0;
            Ogridz[128] = q1;
            Ogridx[129] = -p1;
            Ogridy[129] = 0.0;
            Ogridz[129] = -q1;

            Ogridx[130] = 0.0;
            Ogridy[130] = p1;
            Ogridz[130] = q1;
            Ogridx[131] = 0.0;
            Ogridy[131] = p1;
            Ogridz[131] = -q1;
            Ogridx[132] = 0.0;
            Ogridy[132] = -p1;
            Ogridz[132] = q1;
            Ogridx[133] = 0.0;
            Ogridy[133] = -p1;
            Ogridz[133] = -q1;

            Ogridx[134] = q1;
            Ogridy[134] = p1;
            Ogridz[134] = 0.0;
            Ogridx[135] = q1;
            Ogridy[135] = -p1;
            Ogridz[135] = 0.0;
            Ogridx[136] = -q1;
            Ogridy[136] = p1;
            Ogridz[136] = 0.0;
            Ogridx[137] = -q1;
            Ogridy[137] = -p1;
            Ogridz[137] = 0.0;

            Ogridx[138] = q1;
            Ogridy[138] = 0.0;
            Ogridz[138] = p1;
            Ogridx[139] = q1;
            Ogridy[139] = 0.0;
            Ogridz[139] = -p1;
            Ogridx[140] = -q1;
            Ogridy[140] = 0.0;
            Ogridz[140] = p1;
            Ogridx[141] = -q1;
            Ogridy[141] = 0.0;
            Ogridz[141] = -p1;

            Ogridx[142] = 0.0;
            Ogridy[142] = q1;
            Ogridz[142] = p1;
            Ogridx[143] = 0.0;
            Ogridy[143] = q1;
            Ogridz[143] = -p1;
            Ogridx[144] = 0.0;
            Ogridy[144] = -q1;
            Ogridz[144] = p1;
            Ogridx[145] = 0.0;
            Ogridy[145] = -q1;
            Ogridz[145] = -p1;

            // set dik
            Ogridx[146 + 0] = dr;
            Ogridy[146 + 0] = du;
            Ogridz[146 + 0] = dw;
            Ogridx[146 + 1] = -dr;
            Ogridy[146 + 1] = du;
            Ogridz[146 + 1] = dw;
            Ogridx[146 + 2] = dr;
            Ogridy[146 + 2] = -du;
            Ogridz[146 + 2] = dw;
            Ogridx[146 + 3] = dr;
            Ogridy[146 + 3] = du;
            Ogridz[146 + 3] = -dw;
            Ogridx[146 + 4] = -dr;
            Ogridy[146 + 4] = -du;
            Ogridz[146 + 4] = dw;
            Ogridx[146 + 5] = -dr;
            Ogridy[146 + 5] = du;
            Ogridz[146 + 5] = -dw;
            Ogridx[146 + 6] = dr;
            Ogridy[146 + 6] = -du;
            Ogridz[146 + 6] = -dw;
            Ogridx[146 + 7] = -dr;
            Ogridy[146 + 7] = -du;
            Ogridz[146 + 7] = -dw;

            for (i = 0; i < 8; ++i) {
                Ogridx[146 + 8 + i] = Ogridx[146 + i];
                Ogridy[146 + 8 + i] = Ogridz[146 + i];
                Ogridz[146 + 8 + i] = Ogridy[146 + i];

                Ogridx[146 + 16 + i] = Ogridy[146 + i];
                Ogridy[146 + 16 + i] = Ogridx[146 + i];
                Ogridz[146 + 16 + i] = Ogridz[146 + i];

                Ogridx[146 + 24 + i] = Ogridy[146 + i];
                Ogridy[146 + 24 + i] = Ogridz[146 + i];
                Ogridz[146 + 24 + i] = Ogridx[146 + i];

                Ogridx[146 + 32 + i] = Ogridz[146 + i];
                Ogridy[146 + 32 + i] = Ogridy[146 + i];
                Ogridz[146 + 32 + i] = Ogridx[146 + i];

                Ogridx[146 + 40 + i] = Ogridz[146 + i];
                Ogridy[146 + 40 + i] = Ogridx[146 + i];
                Ogridz[146 + 40 + i] = Ogridy[146 + i];
            }

            // FOR CHECK
            //      for ( i = 0; i < nOgrid; ++i ) {
            //              cout << i  <<  " " <<  Ogridx[i] << "  " <<
            //              Ogridy[i] << " " << Ogridz[i] << "  "
            //                   << (
            //                   Ogridx[i]*Ogridx[i]+Ogridy[i]*Ogridy[i]+Ogridz[i]*Ogridz[i]
            //                   ) << "  " << w[i] << endl;
            //      }
            // END
            // Scale and Shift
            for (i = 0; i < nOgrid; ++i) {
                Ogridx[i] = r0 * Ogridx[i] + corex;
                Ogridy[i] = r0 * Ogridy[i] + corey;
                Ogridz[i] = r0 * Ogridz[i] + corez;
            }

        }  // end-of grid-if

    } else if (grid == 0 || grid == 1) {
        A1 = 0.000599631368862;
        A2 = 0.00737299971862;
        A3 = 0.00721051536014;
        B1 = 0.00757439415905;
        B2 = 0.00675382948631;
        B3 = 0.00711635549312;
        D1 = 0.00699108735330;

        dr = 0.882270011260;
        du = 0.140355381171;
        dw = 0.449332832327;

        CORRP = 1.00;

        bm[0] = 0.974888643677;
        bm[1] = 0.807089818360;
        bm[2] = 0.291298882210;

        // Set Lebedev weight
        // Note CORRP is a temporal correction.
        for (i = 0; i < 6; ++i)
            //                      w[i] *= A1 * 146.0;
            w[i] *= A1 * 146.0 / CORRP;
        for (i = 6; i < 18; ++i)
            //                      w[i] *= A2 * 146.0;
            w[i] *= A2 * 146.0 / CORRP;
        for (i = 18; i < 26; ++i)
            //                      w[i] *= A3 * 146.0;
            w[i] *= A3 * 146.0 / CORRP;
        for (i = 26; i < 50; ++i)
            //                      w[i] *= B1 * 146.0;
            w[i] *= B1 * 146.0 / CORRP;
        for (i = 50; i < 74; ++i)
            //                      w[i] *= B2 * 146.0;
            w[i] *= B2 * 146.0 / CORRP;
        for (i = 74; i < 98; ++i)
            //                      w[i] *= B3 * 146.0;
            w[i] *= B3 * 146.0 / CORRP;
        for (i = 98; i < 146; ++i)
            //                      w[i] *= D1 * 146.0;
            w[i] *= D1 * 146.0 / CORRP;

        // Generate Omega Grid
        // set aik
        Ogridx[0] = 1.0;
        Ogridy[0] = 0.0;
        Ogridz[0] = 0.0;
        Ogridx[1] = 0.0;
        Ogridy[1] = 1.0;
        Ogridz[1] = 0.0;
        Ogridx[2] = 0.0;
        Ogridy[2] = 0.0;
        Ogridz[2] = 1.0;
        Ogridx[3] = -1.0;
        Ogridy[3] = 0.0;
        Ogridz[3] = 0.0;
        Ogridx[4] = 0.0;
        Ogridy[4] = -1.0;
        Ogridz[4] = 0.0;
        Ogridx[5] = 0.0;
        Ogridy[5] = 0.0;
        Ogridz[5] = -1.0;

        Ogridx[6] = SQ1_2;
        Ogridy[6] = SQ1_2;
        Ogridz[6] = 0.0;
        Ogridx[7] = -SQ1_2;
        Ogridy[7] = SQ1_2;
        Ogridz[7] = 0.0;
        Ogridx[8] = SQ1_2;
        Ogridy[8] = -SQ1_2;
        Ogridz[8] = 0.0;
        Ogridx[9] = -SQ1_2;
        Ogridy[9] = -SQ1_2;
        Ogridz[9] = 0.0;

        for (i = 0; i < 4; ++i) {
            Ogridx[6 + 4 + i] = Ogridy[6 + i];
            Ogridy[6 + 4 + i] = Ogridz[6 + i];
            Ogridz[6 + 4 + i] = Ogridx[6 + i];

            Ogridx[6 + 8 + i] = Ogridz[6 + i];
            Ogridy[6 + 8 + i] = Ogridx[6 + i];
            Ogridz[6 + 8 + i] = Ogridy[6 + i];
        }

        Ogridx[18] = SQ1_3;
        Ogridy[18] = SQ1_3;
        Ogridz[18] = SQ1_3;
        Ogridx[19] = -SQ1_3;
        Ogridy[19] = SQ1_3;
        Ogridz[19] = SQ1_3;
        Ogridx[20] = SQ1_3;
        Ogridy[20] = -SQ1_3;
        Ogridz[20] = SQ1_3;
        Ogridx[21] = SQ1_3;
        Ogridy[21] = SQ1_3;
        Ogridz[21] = -SQ1_3;
        Ogridx[22] = -SQ1_3;
        Ogridy[22] = -SQ1_3;
        Ogridz[22] = SQ1_3;
        Ogridx[23] = -SQ1_3;
        Ogridy[23] = SQ1_3;
        Ogridz[23] = -SQ1_3;
        Ogridx[24] = SQ1_3;
        Ogridy[24] = -SQ1_3;
        Ogridz[24] = -SQ1_3;
        Ogridx[25] = -SQ1_3;
        Ogridy[25] = -SQ1_3;
        Ogridz[25] = -SQ1_3;

        // set bik; no need to set cik
        for (i = 0; i < 3; ++i) {
            bl[i] = sqrt(1.0 - bm[i] * bm[i]) / SQ2;

            Ogridx[26 + 24 * i + 0] = bl[i];
            Ogridy[26 + 24 * i + 0] = bl[i];
            Ogridz[26 + 24 * i + 0] = bm[i];
            Ogridx[26 + 24 * i + 1] = -bl[i];
            Ogridy[26 + 24 * i + 1] = bl[i];
            Ogridz[26 + 24 * i + 1] = bm[i];
            Ogridx[26 + 24 * i + 2] = bl[i];
            Ogridy[26 + 24 * i + 2] = -bl[i];
            Ogridz[26 + 24 * i + 2] = bm[i];
            Ogridx[26 + 24 * i + 3] = bl[i];
            Ogridy[26 + 24 * i + 3] = bl[i];
            Ogridz[26 + 24 * i + 3] = -bm[i];
            Ogridx[26 + 24 * i + 4] = -bl[i];
            Ogridy[26 + 24 * i + 4] = -bl[i];
            Ogridz[26 + 24 * i + 4] = bm[i];
            Ogridx[26 + 24 * i + 5] = -bl[i];
            Ogridy[26 + 24 * i + 5] = bl[i];
            Ogridz[26 + 24 * i + 5] = -bm[i];
            Ogridx[26 + 24 * i + 6] = bl[i];
            Ogridy[26 + 24 * i + 6] = -bl[i];
            Ogridz[26 + 24 * i + 6] = -bm[i];
            Ogridx[26 + 24 * i + 7] = -bl[i];
            Ogridy[26 + 24 * i + 7] = -bl[i];
            Ogridz[26 + 24 * i + 7] = -bm[i];
        }

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 8; ++j) {
                Ogridx[26 + 24 * i + 8 + j] = Ogridy[26 + 24 * i + j];
                Ogridy[26 + 24 * i + 8 + j] = Ogridz[26 + 24 * i + j];
                Ogridz[26 + 24 * i + 8 + j] = Ogridx[26 + 24 * i + j];

                Ogridx[26 + 24 * i + 16 + j] = Ogridz[26 + 24 * i + j];
                Ogridy[26 + 24 * i + 16 + j] = Ogridx[26 + 24 * i + j];
                Ogridz[26 + 24 * i + 16 + j] = Ogridy[26 + 24 * i + j];
            }
        }

        // set dik
        Ogridx[98 + 0] = dr;
        Ogridy[98 + 0] = du;
        Ogridz[98 + 0] = dw;
        Ogridx[98 + 1] = -dr;
        Ogridy[98 + 1] = du;
        Ogridz[98 + 1] = dw;
        Ogridx[98 + 2] = dr;
        Ogridy[98 + 2] = -du;
        Ogridz[98 + 2] = dw;
        Ogridx[98 + 3] = dr;
        Ogridy[98 + 3] = du;
        Ogridz[98 + 3] = -dw;
        Ogridx[98 + 4] = -dr;
        Ogridy[98 + 4] = -du;
        Ogridz[98 + 4] = dw;
        Ogridx[98 + 5] = -dr;
        Ogridy[98 + 5] = du;
        Ogridz[98 + 5] = -dw;
        Ogridx[98 + 6] = dr;
        Ogridy[98 + 6] = -du;
        Ogridz[98 + 6] = -dw;
        Ogridx[98 + 7] = -dr;
        Ogridy[98 + 7] = -du;
        Ogridz[98 + 7] = -dw;

        for (i = 0; i < 8; ++i) {
            Ogridx[98 + 8 + i] = Ogridx[98 + i];
            Ogridy[98 + 8 + i] = Ogridz[98 + i];
            Ogridz[98 + 8 + i] = Ogridy[98 + i];

            Ogridx[98 + 16 + i] = Ogridy[98 + i];
            Ogridy[98 + 16 + i] = Ogridx[98 + i];
            Ogridz[98 + 16 + i] = Ogridz[98 + i];

            Ogridx[98 + 24 + i] = Ogridy[98 + i];
            Ogridy[98 + 24 + i] = Ogridz[98 + i];
            Ogridz[98 + 24 + i] = Ogridx[98 + i];

            Ogridx[98 + 32 + i] = Ogridz[98 + i];
            Ogridy[98 + 32 + i] = Ogridy[98 + i];
            Ogridz[98 + 32 + i] = Ogridx[98 + i];

            Ogridx[98 + 40 + i] = Ogridz[98 + i];
            Ogridy[98 + 40 + i] = Ogridx[98 + i];
            Ogridz[98 + 40 + i] = Ogridy[98 + i];
        }

        // FOR CHEOACK
        //      for ( i = 0; i < 146; ++i ) {
        //              cout << " " <<  Ogridx[i] << "  " << Ogridy[i] << "  "
        //              <<
        //              Ogridz[i] << "  "     \
//                   << (
        //                   Ogridx[i]*Ogridx[i]+Ogridy[i]*Ogridy[i]+Ogridz[i]*Ogridz[i]
        //                   ) << endl;
        //      }
        // END
        // Scale and Shift
        for (i = 0; i < nOgrid; ++i) {
            Ogridx[i] = r0 * Ogridx[i] + corex;
            Ogridy[i] = r0 * Ogridy[i] + corey;
            Ogridz[i] = r0 * Ogridz[i] + corez;
        }
    }  // end-of grid-if
}

/*
void points(int grid, int nOgrid, double r0, double* Ogridx, double* Ogridy,
double* Ogridz, double* w, double weight)
{
// Set the grid coordinates in "Ogridr" and the Lebedev weight in "w"
// Get atom core coordinates

        double  corex, corey, corez;
        int     i, j;

        double  A1    = 0.000599631368862;
        double  A2    = 0.00737299971862;
        double  A3    = 0.00721051536014;
        double  B1    = 0.00757439415905;
        double  B2    = 0.00675382948631;
        double  B3    = 0.00711635549312;
        double  D1    = 0.00699108735330;

        double  bm[3];
        double  bl[3];
        double  dr    = 0.882270011260;
        double  du    = 0.140355381171;
        double  dw    = 0.449332832327;

        double  CORRP = 1.00;

        bm[0]         = 0.974888643677;
        bm[1]         = 0.807089818360;
        bm[2]         = 0.291298882210;

        corex = Ogridx[0];
        corey = Ogridy[0];
        corez = Ogridz[0];

        for ( i = 0; i < nOgrid; ++i ) {
                w[i] = weight;
                Ogridx[i] = Ogridy[i] = Ogridz[i] = 0.0;
        }

        if ( grid == -1 ) {

                cout << "Sorry, do not support now." << endl;

        } else if ( grid == 0 || grid == 1 ) {

// Set Lebedev weight
// Note CORRP is a temporal correction.
                for ( i =  0; i <   6; ++i )
//                      w[i] *= A1 * 146.0;
                        w[i] *= A1 * 146.0 / CORRP;
                for ( i =  6; i <  18; ++i )
//                      w[i] *= A2 * 146.0;
                        w[i] *= A2 * 146.0 / CORRP;
                for ( i = 18; i <  26; ++i )
//                      w[i] *= A3 * 146.0;
                        w[i] *= A3 * 146.0 / CORRP;
                for ( i = 26; i <  50; ++i )
//                      w[i] *= B1 * 146.0;
                        w[i] *= B1 * 146.0 / CORRP;
                for ( i = 50; i <  74; ++i )
//                      w[i] *= B2 * 146.0;
                        w[i] *= B2 * 146.0 / CORRP;
                for ( i = 74; i <  98; ++i )
//                      w[i] *= B3 * 146.0;
                        w[i] *= B3 * 146.0 / CORRP;
                for ( i = 98; i < 146; ++i )
//                      w[i] *= D1 * 146.0;
                        w[i] *= D1 * 146.0 / CORRP;

// Generate Omega Grid
                // set aik
                Ogridx[ 0 ] =  1.0;
                Ogridy[ 0 ] =  0.0;
                Ogridz[ 0 ] =  0.0;
                Ogridx[ 1 ] =  0.0;
                Ogridy[ 1 ] =  1.0;
                Ogridz[ 1 ] =  0.0;
                Ogridx[ 2 ] =  0.0;
                Ogridy[ 2 ] =  0.0;
                Ogridz[ 2 ] =  1.0;
                Ogridx[ 3 ] = -1.0;
                Ogridy[ 3 ] =  0.0;
                Ogridz[ 3 ] =  0.0;
                Ogridx[ 4 ] =  0.0;
                Ogridy[ 4 ] = -1.0;
                Ogridz[ 4 ] =  0.0;
                Ogridx[ 5 ] =  0.0;
                Ogridy[ 5 ] =  0.0;
                Ogridz[ 5 ] = -1.0;

                Ogridx[ 6 ] =  SQ1_2;
                Ogridy[ 6 ] =  SQ1_2;
                Ogridz[ 6 ] =  0.0;
                Ogridx[ 7 ] = -SQ1_2;
                Ogridy[ 7 ] =  SQ1_2;
                Ogridz[ 7 ] =  0.0;
                Ogridx[ 8 ] =  SQ1_2;
                Ogridy[ 8 ] = -SQ1_2;
                Ogridz[ 8 ] =  0.0;
                Ogridx[ 9 ] = -SQ1_2;
                Ogridy[ 9 ] = -SQ1_2;
                Ogridz[ 9 ] =  0.0;

                for ( i = 0; i < 4; ++i ) {
                        Ogridx[ 6 + 4 + i ] = Ogridy[ 6 + i ];
                        Ogridy[ 6 + 4 + i ] = Ogridz[ 6 + i ];
                        Ogridz[ 6 + 4 + i ] = Ogridx[ 6 + i ];

                        Ogridx[ 6 + 8 + i ] = Ogridz[ 6 + i ];
                        Ogridy[ 6 + 8 + i ] = Ogridx[ 6 + i ];
                        Ogridz[ 6 + 8 + i ] = Ogridy[ 6 + i ];
                }

                Ogridx[ 18 ] =  SQ1_3;
                Ogridy[ 18 ] =  SQ1_3;
                Ogridz[ 18 ] =  SQ1_3;
                Ogridx[ 19 ] = -SQ1_3;
                Ogridy[ 19 ] =  SQ1_3;
                Ogridz[ 19 ] =  SQ1_3;
                Ogridx[ 20 ] =  SQ1_3;
                Ogridy[ 20 ] = -SQ1_3;
                Ogridz[ 20 ] =  SQ1_3;
                Ogridx[ 21 ] =  SQ1_3;
                Ogridy[ 21 ] =  SQ1_3;
                Ogridz[ 21 ] = -SQ1_3;
                Ogridx[ 22 ] = -SQ1_3;
                Ogridy[ 22 ] = -SQ1_3;
                Ogridz[ 22 ] =  SQ1_3;
                Ogridx[ 23 ] = -SQ1_3;
                Ogridy[ 23 ] =  SQ1_3;
                Ogridz[ 23 ] = -SQ1_3;
                Ogridx[ 24 ] =  SQ1_3;
                Ogridy[ 24 ] = -SQ1_3;
                Ogridz[ 24 ] = -SQ1_3;
                Ogridx[ 25 ] = -SQ1_3;
                Ogridy[ 25 ] = -SQ1_3;
                Ogridz[ 25 ] = -SQ1_3;

                // set bik; no need to set cik
                for ( i = 0; i < 3; ++i ) {
                        bl[i] = sqrt( 1.0 - bm[i] * bm[i] ) / SQ2;

                        Ogridx[ 26 + 24*i + 0 ] =  bl[i];
                        Ogridy[ 26 + 24*i + 0 ] =  bl[i];
                        Ogridz[ 26 + 24*i + 0 ] =  bm[i];
                        Ogridx[ 26 + 24*i + 1 ] = -bl[i];
                        Ogridy[ 26 + 24*i + 1 ] =  bl[i];
                        Ogridz[ 26 + 24*i + 1 ] =  bm[i];
                        Ogridx[ 26 + 24*i + 2 ] =  bl[i];
                        Ogridy[ 26 + 24*i + 2 ] = -bl[i];
                        Ogridz[ 26 + 24*i + 2 ] =  bm[i];
                        Ogridx[ 26 + 24*i + 3 ] =  bl[i];
                        Ogridy[ 26 + 24*i + 3 ] =  bl[i];
                        Ogridz[ 26 + 24*i + 3 ] = -bm[i];
                        Ogridx[ 26 + 24*i + 4 ] = -bl[i];
                        Ogridy[ 26 + 24*i + 4 ] = -bl[i];
                        Ogridz[ 26 + 24*i + 4 ] =  bm[i];
                        Ogridx[ 26 + 24*i + 5 ] = -bl[i];
                        Ogridy[ 26 + 24*i + 5 ] =  bl[i];
                        Ogridz[ 26 + 24*i + 5 ] = -bm[i];
                        Ogridx[ 26 + 24*i + 6 ] =  bl[i];
                        Ogridy[ 26 + 24*i + 6 ] = -bl[i];
                        Ogridz[ 26 + 24*i + 6 ] = -bm[i];
                        Ogridx[ 26 + 24*i + 7 ] = -bl[i];
                        Ogridy[ 26 + 24*i + 7 ] = -bl[i];
                        Ogridz[ 26 + 24*i + 7 ] = -bm[i];
                }

                for ( i = 0; i < 3; ++i ) {
                        for ( j = 0; j < 8; ++j ) {
                                Ogridx[ 26 + 24*i +  8 + j ] = Ogridy[ 26 + 24*i
+ j ]; Ogridy[ 26 + 24*i +  8 + j ] = Ogridz[ 26 + 24*i + j ]; Ogridz[ 26 + 24*i
+  8 + j ] = Ogridx[ 26 + 24*i + j ];

                                Ogridx[ 26 + 24*i + 16 + j ] = Ogridz[ 26 + 24*i
+ j ]; Ogridy[ 26 + 24*i + 16 + j ] = Ogridx[ 26 + 24*i + j ]; Ogridz[ 26 + 24*i
+ 16 + j ] = Ogridy[ 26 + 24*i + j ];
                        }
                }

                // set dik
                Ogridx[ 98 + 0 ] =  dr;
                Ogridy[ 98 + 0 ] =  du;
                Ogridz[ 98 + 0 ] =  dw;
                Ogridx[ 98 + 1 ] = -dr;
                Ogridy[ 98 + 1 ] =  du;
                Ogridz[ 98 + 1 ] =  dw;
                Ogridx[ 98 + 2 ] =  dr;
                Ogridy[ 98 + 2 ] = -du;
                Ogridz[ 98 + 2 ] =  dw;
                Ogridx[ 98 + 3 ] =  dr;
                Ogridy[ 98 + 3 ] =  du;
                Ogridz[ 98 + 3 ] = -dw;
                Ogridx[ 98 + 4 ] = -dr;
                Ogridy[ 98 + 4 ] = -du;
                Ogridz[ 98 + 4 ] =  dw;
                Ogridx[ 98 + 5 ] = -dr;
                Ogridy[ 98 + 5 ] =  du;
                Ogridz[ 98 + 5 ] = -dw;
                Ogridx[ 98 + 6 ] =  dr;
                Ogridy[ 98 + 6 ] = -du;
                Ogridz[ 98 + 6 ] = -dw;
                Ogridx[ 98 + 7 ] = -dr;
                Ogridy[ 98 + 7 ] = -du;
                Ogridz[ 98 + 7 ] = -dw;

                for ( i = 0; i < 8; ++i ) {
                        Ogridx[ 98 +  8 + i ] = Ogridx[ 98 + i ];
                        Ogridy[ 98 +  8 + i ] = Ogridz[ 98 + i ];
                        Ogridz[ 98 +  8 + i ] = Ogridy[ 98 + i ];

                        Ogridx[ 98 + 16 + i ] = Ogridy[ 98 + i ];
                        Ogridy[ 98 + 16 + i ] = Ogridx[ 98 + i ];
                        Ogridz[ 98 + 16 + i ] = Ogridz[ 98 + i ];

                        Ogridx[ 98 + 24 + i ] = Ogridy[ 98 + i ];
                        Ogridy[ 98 + 24 + i ] = Ogridz[ 98 + i ];
                        Ogridz[ 98 + 24 + i ] = Ogridx[ 98 + i ];

                        Ogridx[ 98 + 32 + i ] = Ogridz[ 98 + i ];
                        Ogridy[ 98 + 32 + i ] = Ogridy[ 98 + i ];
                        Ogridz[ 98 + 32 + i ] = Ogridx[ 98 + i ];

                        Ogridx[ 98 + 40 + i ] = Ogridz[ 98 + i ];
                        Ogridy[ 98 + 40 + i ] = Ogridx[ 98 + i ];
                        Ogridz[ 98 + 40 + i ] = Ogridy[ 98 + i ];
                }

// FOR CHECK
//      for ( i = 0; i < 146; ++i ) {
//              cout << " " <<  Ogridx[i] << "  " << Ogridy[i] << "  " <<
Ogridz[i] << "  "     \
//                   << (
Ogridx[i]*Ogridx[i]+Ogridy[i]*Ogridy[i]+Ogridz[i]*Ogridz[i] ) << endl;
//      }
// END
// Scale and Shift
                for ( i = 0; i < nOgrid; ++i ) {
                        Ogridx[i] = r0 * Ogridx[i] + corex;
                        Ogridy[i] = r0 * Ogridy[i] + corey;
                        Ogridz[i] = r0 * Ogridz[i] + corez;
                }

        } else if ( grid == 2 ) {

                cout << "Sorry, do not support now." << endl;

        }                                       // end-of grid-if

}
*/

double polfunc(double z) {
    // calculation polarization function

    double f, p, m;

    p = pow(1.0 + z, F43);
    m = pow(1.0 - z, F43);
    f = 0.5 * (p + m - 2.0) / (R3_2 - 1.0);

    return f;
}

double poldrfunc(double z) {
    // calculation drived polarization function

    double f, p, m;

    p = pow(1.0 + z, F13);
    m = pow(1.0 - z, F13);
    f = 2.0 * F13 * (p - m) / (R3_2 - 1.0);

    return f;
}

double HLfunc(double z) {
    // calculation Hedin & Lundqvist function
    // it is also used on JMW and GL.

    double fz;

    fz = (1.0 + z * z * z) * log(1.0 + 1.0 / z);
    fz += z / 2.0 - z * z - F13;

    return fz;
}

double HLdrfunc(double z) {
    // calculation Hedin & Lundqvist drived function
    // it is also used on JMW and GL.

    double dfz;

    dfz = 3.0 * z * z * log(1.0 + 1.0 / z);
    dfz += 1.5 - 3.0 * z - 1.0 / z;

    return dfz;
}

double VWNPfunc(double z) {
    // calculation VWN P function

    double fz, x, XPx, atQ;
    double aP = -0.10498;
    double bP = 3.72744;
    double CP = 12.9352;
    double QP = 6.15199081976;   // QP  = sqrt( 4.0*CP - bP*bP )
    double XPa = 12.5549141492;  // XPa = aP*aP + bP*aP + CP
                                 //      double  aP  = -0.409286;
                                 //      double  bP  = 13.0720;
                                 //      double  CP  = 42.7198;
                                 //      double  QP, XPa;
                                 //      QP  = sqrt( 4.0*CP - bP*bP );
                                 //      XPa = aP*aP + bP*aP + CP;

    x = sqrt(z);
    XPx = x * x + bP * x + CP;
    atQ = atan(QP / (2.0 * x + bP));

    fz = log((x - aP) * (x - aP) / XPx);
    fz += 2.0 * (bP + 2.0 * aP) * atQ / QP;
    fz *= -bP * aP / XPa;
    fz += log(x * x / XPx) + 2.0 * bP * atQ / QP;

    return fz;
}

double VWNPdrfunc(double z) {
    // calculation VWN P derived function

    double dfz, x, XPx, xbP, xQQ, xxX;
    double aP = -0.10498;
    double bP = 3.72744;
    double CP = 12.9352;
    double QP = 6.15199081976;   // QP  = sqrt( 4.0*CP - bP*bP )
    double XPa = 12.5549141492;  // XPa = aP*aP + bP*aP + CP
                                 //      double  aP  = -0.409286;
                                 //      double  bP  = 13.0720;
                                 //      double  CP  = 42.7198;
                                 //      double  QP, XPa;
                                 //      QP  = sqrt( 4.0*CP - bP*bP );
                                 //      XPa = aP*aP + bP*aP + CP;

    x = sqrt(z);
    XPx = x * x + bP * x + CP;
    xbP = 2.0 * x + bP;
    xQQ = 4.0 * x / (QP * QP + xbP * xbP);
    xxX = xbP * x / XPx;

    dfz = xQQ * (bP + 2.0 * aP);
    dfz -= 2.0 * x / (x - aP) + xxX;
    dfz *= bP * aP / XPa;
    dfz += 2.0 - xxX - xQQ * bP;

    return dfz;
}

double VWNFfunc(double z) {
    // calculation VWN F function

    double fz, x, XPx, atQ;
    double aP = -0.32500;
    double bP = 7.06042;
    double CP = 18.0578;
    double QP = 4.73092690956;  // QP  = sqrt( 4.0*CP - bP*bP )
    double XPa = 15.8687885;    // XPa = aP*aP + bP*aP + CP
                                //      double  aP  =  -0.743294;
                                //      double  bP  =  20.1231;
                                //      double  CP  = 101.578;
                                //      double  QP, XPa;
                                //      QP  = sqrt( 4.0*CP - bP*bP );
                                //      XPa = aP*aP + bP*aP + CP;

    x = sqrt(z);
    XPx = x * x + bP * x + CP;
    atQ = atan(QP / (2.0 * x + bP)) / QP;

    fz = log((x - aP) * (x - aP) / XPx);
    fz += 2.0 * (bP + 2.0 * aP) * atQ;
    fz *= -bP * aP / XPa;
    fz += log(x * x / XPx) + 2.0 * bP * atQ;

    return fz;
}

double VWNFdrfunc(double z) {
    // calculation VWN F derived function

    double dfz, x, XPx, xbP, xQQ, xxX;
    double aP = -0.32500;
    double bP = 7.06042;
    double CP = 18.0578;
    double QP = 4.73092690956;  // QP  = sqrt( 4.0*CP - bP*bP )
    double XPa = 15.8687885;    // XPa = aP*aP + bP*aP + CP
                                //      double  aP  =  -0.743294;
                                //      double  bP  =  20.1231;
                                //      double  CP  = 101.578;
                                //      double  QP, XPa;
                                //      QP  = sqrt( 4.0*CP - bP*bP );
                                //      XPa = aP*aP + bP*aP + CP;

    x = sqrt(z);
    XPx = x * x + bP * x + CP;
    xbP = 2.0 * x + bP;
    xQQ = 4.0 * x / (QP * QP + xbP * xbP);
    xxX = xbP * x / XPx;

    dfz = xQQ * (bP + 2.0 * aP);
    dfz -= 2.0 * x / (x - aP) + xxX;
    dfz *= bP * aP / XPa;
    dfz += 2.0 - xxX - xQQ * bP;

    return dfz;
}

double B88func(double RouA, double RouB, double XA, double XB) {
    // calculate b88 functional
    double GCex, GCexA, GCexB;
    double b = 0.0042;
    double c = 0.9305257363491002;  // 3/2 * pow( 3/(4*pai), 1/3 )

    GCexA = -b * XA * XA / (1.0 + 6.0 * b * XA * asinh(XA));
    GCexB = -b * XB * XB / (1.0 + 6.0 * b * XB * asinh(XB));
    GCex = pow(RouA, F43) * (GCexA - c) + pow(RouB, F43) * (GCexB - c);

    return GCex;
}

double DB88func(double RouA, double RouB, double XA, double XB) {
    // calculate b88 functional except const
    double GCex, GCexA, GCexB;
    double b = 0.0042;

    GCexA = -b * XA * XA / (1.0 + 6.0 * b * XA * asinh(XA));
    GCexB = -b * XB * XB / (1.0 + 6.0 * b * XB * asinh(XB));
    GCex = pow(RouA, F43) * GCexA + pow(RouB, F43) * GCexB;

    return GCex;
}

double B88dfunc(double Rou, double X, double G, double gRx, double gRy,
                double gRz, double g, double ggx, double ggy, double ggz) {
    // calculate first derivative of b88 functional
    double dGCex;
    double gX, dgX;
    double dB_dR, dB_dG;
    double XX, ish, gRgg, osbXi;
    double b = 0.0042;
    double bb = 0.00001764;         // b * b
    double c = 0.9305257363491002;  // 3/2 * pow( 3/(4*pai), 1/3 )

    XX = X * X;
    ish = asinh(X);
    gRgg = gRx * ggx + gRy * ggy + gRz * ggz;
    osbXi = 1.0 + 6.0 * b * X * ish;

    gX = -b * XX / osbXi - c;
    dgX = (6.0 * bb * XX * (X / sqrt(XX + 1.0) - ish) - 2.0 * b * X) /
          (osbXi * osbXi);

    dB_dR = F43 * pow(Rou, F13) * (gX - X * dgX);
    dB_dG = 0.5 * dgX / sqrt(G);

    dGCex = dB_dR * g + 2.0 * dB_dG * gRgg;

    return dGCex;
}

double DB88dfunc(double Rou, double X, double G, double gRx, double gRy,
                 double gRz, double g, double ggx, double ggy, double ggz) {
    // calculate first derivative of b88 functional except const
    double dGCex;
    double gX, dgX;
    double dB_dR, dB_dG;
    double XX, ish, gRgg, osbXi;
    double b = 0.0042;
    double bb = 0.00001764;  // b * b

    XX = X * X;
    ish = asinh(X);
    gRgg = gRx * ggx + gRy * ggy + gRz * ggz;
    osbXi = 1.0 + 6.0 * b * X * ish;

    gX = -b * XX / osbXi;
    dgX = (6.0 * bb * XX * (X / sqrt(XX + 1.0) - ish) - 2.0 * b * X) /
          (osbXi * osbXi);

    dB_dR = F43 * pow(Rou, F13) * (gX - X * dgX);
    dB_dG = 0.5 * dgX / sqrt(G);

    dGCex = dB_dR * g + 2.0 * dB_dG * gRgg;

    return dGCex;
}

double G96func(double RouA, double RouB, double XA, double XB) {
    // calculate g96 functional
    double GCex, GCexA, GCexB;
    double B = 0.0072992700729927;  // 1.0/137.0;
    double c = 0.9305257363491002;  // 3/2 * pow( 3/(4*pai), 1/3 )

    GCexA = -B * pow(XA, 1.5);
    GCexB = -B * pow(XB, 1.5);
    GCex = pow(RouA, F43) * (GCexA - c) + pow(RouB, F43) * (GCexB - c);

    return GCex;
}

double DG96func(double RouA, double RouB, double XA, double XB) {
    // calculate g96 functional except const
    double GCex, GCexA, GCexB;
    double B = 0.0072992700729927;  // 1.0/137.0;

    GCexA = -B * pow(XA, 1.5);
    GCexB = -B * pow(XB, 1.5);
    GCex = pow(RouA, F43) * GCexA + pow(RouB, F43) * GCexB;

    return GCex;
}

double G96dfunc(double Rou, double X, double G, double gRx, double gRy,
                double gRz, double g, double ggx, double ggy, double ggz) {
    // calculate first derivative of g96 functional
    double dGCex;
    double gX, dgX;
    double dB_dR, dB_dG;
    double gRgg;
    double B = 0.0072992700729927;  // 1.0/137.0;
    double c = 0.9305257363491002;  // 3/2 * pow( 3/(4*pai), 1/3 )

    gRgg = gRx * ggx + gRy * ggy + gRz * ggz;

    gX = -B * pow(X, 1.5) - c;
    dgX = -B * 1.5 * sqrt(X);

    dB_dR = F43 * pow(Rou, F13) * (gX - X * dgX);
    dB_dG = 0.5 * dgX / sqrt(G);

    dGCex = dB_dR * g + 2.0 * dB_dG * gRgg;

    return dGCex;
}

double DG96dfunc(double Rou, double X, double G, double gRx, double gRy,
                 double gRz, double g, double ggx, double ggy, double ggz) {
    // calculate first derivative of g96 functional except
    double dGCex;
    double gX, dgX;
    double dB_dR, dB_dG;
    double gRgg;
    double B = 0.0072992700729927;  // 1.0/137.0;

    gRgg = gRx * ggx + gRy * ggy + gRz * ggz;

    gX = -B * pow(X, 1.5);
    dgX = -B * 1.5 * sqrt(X);

    dB_dR = F43 * pow(Rou, F13) * (gX - X * dgX);
    dB_dG = 0.5 * dgX / sqrt(G);

    dGCex = dB_dR * g + 2.0 * dB_dG * gRgg;

    return dGCex;
}

double LYPfunc(double TRou, double RouA, double RouB, double GAA, double GAB,
               double GBB) {
    // calculate lyp functional

    double a, b, c, d;
    double Cfp =
        36.4623989787647670;  // pow( 2, 11/3 ) * 3/10 * pow( 3*pai*pai, 2/3 )
    double Omega, Delta;
    double Rou13, abO, RouAB, odR;
    double F19, ttDelta, Delta11;
    double dLYP_dGAA, dLYP_dGAB, dLYP_dGBB;
    double LYP1, LYP2, LYP3;
    double LYP;

    a = 0.04918;  // Case1 of Parameter   :Default
    b = 0.132;
    c = 0.2533;
    d = 0.349;
    //    a = 0.049;                // Case2
    //    b = 0.108;
    //    c = 0.24;
    //    d = 0.342;

    Rou13 = pow(TRou, -F13);
    RouAB = RouA * RouB;
    odR = 1.0 + d * Rou13;
    Omega = exp(-c * Rou13) * pow(TRou, -11.0 / 3.0) / odR;
    Delta = c * Rou13 + d * Rou13 / odR;
    abO = a * b * Omega;

    F19 = 1.0 / 9.0;
    ttDelta = 1.0 - 3.0 * Delta;
    Delta11 = Delta - 11.0;

    dLYP_dGAA =
        -abO * (F19 * RouAB * (ttDelta - Delta11 * RouA / TRou) - RouB * RouB);
    dLYP_dGBB =
        -abO * (F19 * RouAB * (ttDelta - Delta11 * RouB / TRou) - RouA * RouA);
    dLYP_dGAB = -abO * (F19 * RouAB * (47.0 - 7.0 * Delta) - F43 * TRou * TRou);

    LYP1 = -4.0 * a * RouAB / (odR * TRou);
    LYP2 = -Cfp * abO * RouAB * (pow(RouA, 8.0 / 3.0) + pow(RouB, 8.0 / 3.0));
    LYP3 = dLYP_dGAA * GAA + dLYP_dGAB * GAB + dLYP_dGBB * GBB;
    LYP = LYP1 + LYP2 + LYP3;

    return LYP;
}

double LYPdfunc(double TRou, double RouA, double RouB, double GAA, double GAB,
                double GBB, double gRAx, double gRAy, double gRAz, double gRBx,
                double gRBy, double gRBz, double g, double ggx, double ggy,
                double ggz) {
    // calculate first derivative of lyp functional

    double a, b, c, d, F19, F83;
    double Cfp =
        36.4623989787647670;  // pow( 2, 11/3 ) * 3/10 * pow( 3*pai*pai, 2/3 )
    double Omega, Delta;
    double dOmega, dDelta, abdO, dOmegap;
    double ttDelta, Delta11, fssDelta;
    double Rou13, Rou43, abO, RouAA, RouAB, RouBB, odR;
    double RouA83, RouB83, RouA_R, RouB_R;
    double gRAgg, gRBgg;
    double dLYP1, dLYP2, dLYP3;
    double dLYP_dGAA, dLYP_dGAB, dLYP_dGBB;
    double ddLYP_dRdGAA, ddLYP_dRdGAB, ddLYP_dRdGBB;
    double dLYP_dR;
    double dLYP;

    a = 0.04918;  // Case1 of Parameter   :Default
    b = 0.132;
    c = 0.2533;
    d = 0.349;
    //    a = 0.049;                // Case2
    //    b = 0.108;
    //    c = 0.24;
    //    d = 0.342;

    F19 = 1.0 / 9.0;
    F83 = 8.0 / 3.0;
    gRAgg = gRAx * ggx + gRAy * ggy + gRAz * ggz;
    gRBgg = gRBx * ggx + gRBy * ggy + gRBz * ggz;

    Rou13 = pow(TRou, -F13);
    Rou43 = pow(TRou, -F43);
    RouA83 = pow(RouA, F83);
    RouB83 = pow(RouB, F83);
    RouAA = RouA * RouA;
    RouBB = RouB * RouB;
    RouAB = RouA * RouB;
    RouA_R = RouA / TRou;
    RouB_R = RouB / TRou;
    odR = 1.0 + d * Rou13;

    Omega = exp(-c * Rou13) * pow(TRou, -11.0 / 3.0) / odR;
    Delta = c * Rou13 + d * Rou13 / odR;
    abO = a * b * Omega;
    ttDelta = 1.0 - 3.0 * Delta;
    Delta11 = Delta - 11.0;
    fssDelta = 47.0 - 7.0 * Delta;

    dOmegap = -F13 * Rou43 * (11.0 / Rou13 - c - d / odR);
    dOmega = dOmegap * Omega;
    dDelta = F13 * (d * d * pow(TRou, -5.0 / 3.0) / (odR * odR) - Delta / TRou);
    abdO = a * b * dOmega;

    dLYP_dGAA = -abO * (F19 * RouAB * (ttDelta - Delta11 * RouA_R) - RouBB);
    dLYP_dGBB = -abO * (F19 * RouAB * (ttDelta - Delta11 * RouB_R) - RouAA);
    dLYP_dGAB = -abO * (F19 * RouAB * fssDelta - F43 * TRou * TRou);

    ddLYP_dRdGAA =
        dOmegap * dLYP_dGAA -
        abO *
            (F19 * RouB * (ttDelta - Delta11 * RouA_R) -
             F19 * RouAB * ((3.0 + RouA_R) * dDelta + Delta11 * RouB_R / TRou));
    ddLYP_dRdGAB =
        dOmegap * dLYP_dGAB -
        abO * (F19 * RouB * fssDelta - 7.0 * F19 * RouAB * dDelta - F83 * TRou);
    ddLYP_dRdGBB =
        dOmegap * dLYP_dGBB -
        abO *
            (F19 * RouB * (ttDelta - Delta11 * RouB_R) -
             F19 * RouAB * ((3.0 + RouB_R) * dDelta - Delta11 * RouB_R / TRou) -
             2.0 * RouA);

    dLYP1 = -4.0 * a * RouAB *
            (F13 * d * Rou43 / odR + 1.0 / RouA - 1.0 / TRou) / (odR * TRou);
    dLYP2 = -Cfp * abdO * RouAB * (RouA83 + RouB83) -
            Cfp * abO * RouB * (11.0 / 3.0 * RouA83 + RouB83);
    dLYP3 = ddLYP_dRdGAA * GAA + ddLYP_dRdGAB * GAB + ddLYP_dRdGBB * GBB;
    dLYP_dR = dLYP1 + dLYP2 + dLYP3;

    dLYP = dLYP_dR * g + 2.0 * dLYP_dGAA * gRAgg + dLYP_dGAB * gRBgg;

    return dLYP;
}

/* end of file */ /* DfGrid.cxx */
