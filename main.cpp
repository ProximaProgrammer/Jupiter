#include <iostream>
#include "regimes.h"

#include <vector>
#include <array>
#include <tuple>
#include <any>
using namespace std;

/*
smallest scale height is 15-20 km, while radius of Jupiter is 70000 km, so we prefer N ≥ 4500 if using uniform spherical coordinate cell division
with spherical coordinates, cells near poles become very tiny and lead to numerical instability, 
so use averaging (set all 0 < \theta < \theta_min cells to the same values when one of those cells' values changes, and use the total polar region's volume for calculations instead of one cell)
But see this: https://www.google.com/search?aep=11&q=shortest+scale+height+or+characteristic+length+of+any+phenomena+inside+jupiter+%28shortest+characteristic+scale+where+%E2%88%86x%2Fx+%3C%3C+1+does+not+hold+for+some+variable+x%29&mstk=AUtExfD1bxxmO10akJd4MsbrFn0KhjS8aCB6A2INAwed_t4EEB9TZaVu4hNKULvb2uoYrZyECm6uyOs6zWJouo98S7mE-3z5mcGVgPR27R8GX2O_QP4P00ug9kGGslb-k_1xgD4fDT7JVAjJCbjv3IZROAxlXuLqWFtGnPy3kUpGo_eD2w2bDMVGnf3xVTtGJC17zUi3VG6Y3lXK0lIkJRNoKEdw8Fe2qoujH3UUmtgGsSz6v5sTrMdMSOBXRS5YtbjLlCpCz4bJJsAHVQDpfhRAO2gh0Hcu_C62YITXLizdwUo9I95SU3mnGHcYmTpTB-Y_3EwoBSilYhabaJt9y_eMvPz_IqULAtfY9ru7ard3KTsia6WScIrfDkUhvfTSr2DzOY5Jkp-O_5oV_zpeXcbGIcmwEiC2T3OxSQ&csuir=1&mtid=ml07atDGGtfD0PEPx8zCgAE&udm=50

look into Io-Jupiter and cometary interaction where a lot of ions are channeled into Jupiter and contribute to its magnetic field (may be significant)
*/

struct LineParameter {
    string type;            // Example: "Voigt" or "Synchrotron"
    vector<double> values;  // Store your numerical parameters here
};

int main() {
    //boundary conditions: outer layer emit()s to simulate received sunlight
    //initial conditions: estimated conditions of Jupiter
    const int N = 10; 
    array<array<array<region,N>,N>,N> cells;
    //[some cpp vector] cont_spectrum = ....
    vector<LineParameter> disc_spectrum; //vector of lines with their characterizing parameters (type of line should be a string) Example: {{"Voigt",{parameters}}, {"Synchrotron",{parameters}}}

    //note: a copy of the cells array is made to store the next timestamp's data, calculated on the current data, then the 'older' array is destroyed (ENSURE ITS MEMORY IS CLEANED)
    //since light timescales are far shorter than matter change timescales, we use a quasi-static model where RT is computed everywhere (intermediate array?) and THEN the new cells array based on that (light data is simply copied from the intermediate array)
    //check that in the quasi-static model, we can get away with 1 light iteration for each matter iteration, or whether we need numerous sub-iterations of RT 'towards equilibrium' between each main (matter) iteration
}

//remember that we have the google colab and remote SSH extensions installed ;)