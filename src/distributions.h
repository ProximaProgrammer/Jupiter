#include <iostream>
#include <array>
#include <vector>
#include <tuple>
#include <cmath>
using namespace std;
const double sensor_bin_size = 8.92857143e-9; //AVERAGE wavelength bandwidth of each bin (just for now), since apparently bandwidths & centers vary across the bins
constexpr double c = 299792458.0; //most memory efficient, and using double avoids int truncation issues
const int RESOLUTION = 75000; //logarithmic wavelength bins from 10^-11 m to 10^4 m, with 0.002 dex resolution (500 bins per power of 10)

double blackbody(double lambda, array<double,2> other_params) {
    return 10000; //placeholder for testing
    //NEED TO FACTOR IN EMISSIVITY!

    double temperature = other_params[0];
    double emissivity = other_params[1];
    return sensor_bin_size * 1.19104297e-16 * pow(lambda,-5) / (exp(0.0143877688/(lambda*temperature)) - 1); //multiply bin_size to approximate integral of spectral density over each bin
}

enum chemical_species {//still have to add more, must be region type dependent (move computation of spectral lines to regimes.h)
    H2, He, NH3, H2O, H2S, PH3, NH4SH, HCN, //main stuff
    CH4, C2H6, C2H2, C6H6, CO2}; //organics

//each row is a different chemical with its spectral lines in tuples of the form (frequency f_0, characteristic width gamma)
//NOTE: THIS IS ONLY A PLACEHOLDER; CHANGE LATER AND MATCH TO ENUM OF CHEMICALS (count different phases/forms of the same thing as different chemicals if there are different energies?)
array<vector<tuple<double,double>>,3> lines_example = {{
    {{1e14, 1e-20}, {5e14, 3.233e-18}},
    {{0,0}, {0,0}, {0,0}, {0,0}, {0,0}},
    {{0,0}, {0,0}}
}};
//above, we define a 'dictionary' to plug in a chemical as the index (recall enums evaluate to int so it's really an array)

double STORE_lorentz()

vector<double> OUTPUT_lorentz(double lambda, array<double,4> other_params) { //this explicitly converts characterizing parameters to an observable spectrum; should only be used at very last step (instrument detection)
    //check whether you have to apply shifting and broadening on lines before or after lorentz function
    //CONSIDER: say ions are moving at a velocity different from the surrounding medium. Then how do we deal with Shift() and the two distinct velocities \vec{v} within the same cell?


    
    double intensity_0 = other_params[0];
    double chemical = other_params[1]; //use an enum w/ names so it is actually an int, and the int will automatically cast to double in other_params (need type of all params to be double)
    //^^shouldn't we be considering weighted avg of all present chemicals?? why just one
    double number_density = other_params[2];
    array<vector<tuple<double,double>>,3> lines = other_params[3] ? other_params[3] : lines_example;
    vector<double> spectrum;
   
    // FIX THIS!
    // for (tuple<double,double,double> pair : lines[chemical]) {
    //     double f_0 = get<0>(pair);
    //     double gamma = get<1>(pair);
    //     double intensity_per_number_density = get<2>(pair);
    //     for (double lambda: )
    //     spectrum.push_back(intensity_0 * 7593831.67*gamma/((c/lambda-f_0)*(c/lambda-f_0) + 0.0063325739*gamma*gamma)); //distribution d\lambda, not df
    // }
    return vector<double>();
}



//-----------corner of shame------------//

//return intensity_0 * 7593831.67*gamma/((c/lambda-f_0)*(c/lambda-f_0) + 0.0063325739*gamma*gamma); //distribution d\lambda, not df