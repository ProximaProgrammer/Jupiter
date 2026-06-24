// #include <iostream>
// #include <array>
// #include "distributions.h"
// using namespace std;

// //for computational code use parallelization with the RTX (parallelization may be also efficient to do in python to process multiple files at once)

// //in SPE RDR files, intensities are not normalized but rather true and calibrated spectral radiance
// //to fetch individual SPE RDR file bandwiths and bin center wavelengths, search here??: https://atmos.nmsu.edu/PDS/data/PDS4/juno_jiram_bundle/calibration/

// // changing 'Func' to an explicit function pointer type that uses N
// template <size_t N>
// array<double,336> generate_spectrum(double (*func)(double, const std::array<double, N>), const std::array<double, N> other_params) {
//     std::array<double, 336> spectrum;
    
//     for (int i = 0; i < 336; i++) { 
//         spectrum[i] = func(2e-6 + sensor_bin_size * (i * 1.0 + 0.5), other_params); //bin_size is known when we call 'include "distributions.h"'
//     }
//     return spectrum;
// }

// array<double,RESOLUTION> generate_cont_spectrum(double (*func)(double, const std::array<double, N>), const std::array<double, N> other_params) {
//     std::array<double, RESOLUTION> spectrum;
    
//     for (int i = 0; i < RESOLUTION; i++) { 
//         spectrum[i] = func(2e-6 + sensor_bin_size * (i * 1.0 + 0.5), other_params); //bin_size is known when we call 'include "distributions.h"'
//     }
//     return spectrum;
// }

// void ADD(array<double,RESOLUTION> &arr1, array<double,RESOLUTION> arr2) { //makes alias for arr1 to modify original array and add arr2 to it
//     for (int i=0;i<RESOLUTION;i++) {
//         arr1[i] += arr2[i];
//     }
// }

// //each unit cell dV_0 is an instantiated object of type region
// struct region {
//     double P; //pressure
//     double T; //temperature
//     double density;

//     array<double,3> E;
//     array<double,3> B;
//     double refraction; //index of refraction
//     double electron_fraction;
//     double thermal_diffusivity; //alpha = conductivity/(density x heat capacity)      is the constant in the Heat Equation PDE
//     static constexpr int N_species = 3; //placeholder
//     array<double,N_species> composition_by_mass = {}; //placeholder (compiler wants a guarantee that N_species is constant and array length won't change idk)
//     //can compute molecular mass from composition

//     double emissivity[RESOLUTION]; //emissivity = absorptivity for any wavelength
//     double refractive_index[RESOLUTION]; //can't this be derived with simple math from driven oscillator model?
//     array<double,RESOLUTION> final_emission_spectrum = {}; //initialized to zeroes
//     /*
//     we need a continuous spectral density function for intermediate calculations/basically everything; but we can't just add up evaluated spectrum functions indefinitely
//     this is because we can't lose fine details in the resolution but also don't want to accumulate complexity (just add up in each step)
//     note that we only have a limited set of equations of state, so our spectrum is a superposition line features, thermal distribution, and a few more
//     solution: embed components as a set of width/height/etc. and easily alter these when doppler shift and other modifications apply. (ex.) Lorentz distribution --> store gamma, total line intensity, center wavelength
//     A Voigt profile convolutes a lorentzian with a gaussian (doppler and other broadening), so store Voigt profiles instead of lorentzians
//     ------------------------------------------
//     A/(1 + x^2) + B/(1 + x^2) = (A+B)/(1 + x^2)
//     For x << 1, 1/(1 + x^2) + 1/(1 + 4x^2) = 2/(1 + 5/2 x^2). Note that 5/2 = HM(1,4).
//     For B ~< x << 1 ~ A, 1/(1 + x^2) + A/(1 + (x-B)^2) = (1 + A + (A+1)x^2 - 2Bx + B^2)/(1 + (2x^2 - 2Bx + B^2)) = (1+A)/(1 + x^2) to first order; for B >> x we simply have separate lines and don't combine anything
//     */
//     array<double,RESOLUTION> continuous_emission_spectrum = {}; //also initialized to zeroes (may change once boundary conditions are inserted)

//     virtual void emit() {
//         //THESE ARE REFERENCE FUNCTIONS, MEANT TO BE USED SELECTIVELY WITH RADIATIVE TRANSPORT AND PLUGGED WITH PARAMETERS
//         //can neglect ion charge-exchange X-rays

//         //---------SOURCES---------
//         //thermal blackbody distribution             __characterized by__     T, emissivity(lambda)
//         ADD(continuous_emission_spectrum, generate_spectrum(blackbody, {T})); //NOTE: emissivity E(lambda)=1 cannot be assumed
//         //Synchrotron Radiation                      __characterized by__     T_e, \vec{B}, index p of power law for relativistic electron energy distribution, number density of emitting electrons
//         //Collision-Induced Emission                 __characterized by__     T, colliding species density, composition (determines induced dipole moment as a function of intermolecular distance)
//         //H3+ (trihydrogion cation) emission         __characterized by__     T_rotational, T_vibrational, column density, ionization rate
//         ///--> SEE: https://tinyurl.com/look-at-this-muffins and https://www.exomol.com/data/molecules/H3_p/1H3_p/NMT/ (≈continuous spectrum)
//         //Jovian DAM & HOM radio burst               __characterized by__     electron cyclotron frequency, plasma frequency, pitch angle distribution of energetic electrons
//         //FUV/EUV electron impact excitation of H/H2 __characterized by__     electron energy spectrum, number densities of H/H2, quantum excitation cross sections

//         //Chemical emission (ordinary) LINES         __characterized by__     T, number density
//         //Auroral LINES                              __characterized by__     incident particle energy, neutral atmospheric density, composition, excitation cross-sections
//         /// Io and comets contribute ions, and there is proton plasma too. The emission is anisotropic (due to heavy doppler shifting). Secondary electron impact lines are included in FUV/EUV electron impact emission (physically indistinguishable); see above
//         ///--> SEE: https://www.google.com/search?client=firefox-b-1-d&q=just+as+how+T+and+emissivity%28%5Clambda%29+determine+a+blackbody+distribution%2C+which+parameters+determine+each+of+the+following+spectra%3A+%2F%2FSynchrotron+Radiation+++++++++%2F%2FAuroral+and+Atmospheric+emission+++++++++%2F%2FCollision-Induced+Emission+++++++++%2F%2FH3%2B+%28trihydrogion+cation%29+line+emissions+++++++++%2F%2FJovian+Decametric+%28DAM%29+and+Hectometric+%28HOM%29+radio+burst+++++++++%2F%2FFUV%2FEUV+electron+impact+escitation+of+H%2FH2&udm=50&fbs=ABfTbFVyMZGZf1hfvX9uKjN_-G8c4u0nXx4bEIpwm1lnNH832a9BVCEiB2iPJNekNderQwJGZIG7YID1eBGNWasq2rzBJP_8s1JrXgV7_edZdNiJivgiwJ-uqwgGt4oPbUYp3eFhFdsvSh2uz9KnDVadc_csa1enu2oNy7Ig_i-MHRvNmxvUz5GkAj3ypS8MlqL-FjUf5qSfKVDsQPc_Bif9TfzIECGctw&aep=10&ntc=1&mstk=AUtExfDVfBJEP_XqU26-c7S1SlHOAMrl9wYy0ZSYxyHnC2NyY6AVzNVF46ZqyufVd6LMrWFOjiyauFJhKCvHRuToWJeQgN87PGV2dXGSJD6c-uOxT9E1bBqZUcwhgFf54l41tXaNNA7KrP1GBMur-QN0s2H5k2d4AV2NNThX8ZlSwH6exmstLEWlfAmA4m5xSDniv8vGoiyaallHTDnVwWtnIP89zw8PfZb8CXIsPTDMY51ZtPpj18wjrHtv3LbYRjysa_In1y5qzsdGYosJEtwSZxnocZ9zOYb_or2xHw2ynjySEDtjwDOnPjTaBvt-ac9tVm0mYedr1GGkRzltxS4wGtbbWR_dJ34Hcw&aioh=3&csuir=1&mtid=FDQ5auWEFI-xur8Pj9a7mAg

//         //---------MODIFIERS---------
//         //natural broadening
//         //nest modifiers. (ex.) only doppler shift: ADD(continuous_emission_spectrum, doppler_shift(a, b, c, generate_spectrum(lorentz, {T, the other parameters required})));
//         //doppler broadening - normal, bulk, and rotational
//         //pressure broadening (consists of several factors)
//         //zeeman splitting (significant in magnetically active regions such as core)
//     }
//     virtual void interact() { //be wary of polarized radiation
//         //absorption and scattering (transmission & reflection)
//         //absorption ___
//         //For RT, see this: https://chatgpt.com/share/6a36cd43-5e28-83ed-aecb-68cc8550552a
        
//         /*
//         most general equation of radiative transfer: https://en.wikipedia.org/wiki/Radiative_transfer#The_equation_of_radiative_transfer
//         we need to define I(f, \vec{r}, propagation direction \vec{k}) cellwise and plug it into Radiative Transfer models
//         the composition impacts opacity and optical properties, which varies between region types but also WITHIN a region type (we may not be able to take a simple average)
//         CONSIDER: would monte carlo on photons (w/ random wavevector or something) work better than cellwise Radiative Transfer?
//         */
//         //we need interaction to be wavelength-dependent, still only considering the 336 discrete bins, since that's all JIRAM observes
//     }
//     virtual void change() {
//         //NOTE: these should all be inherited by region types, with equations defined as functions outside any struct
//         //Continuity equation(s)
//         //Navier-Stokes (fluid dynamics, modify with electromagnetic terms as needed)
//         //all other general equations of state; though we require a different set of equations of state depending on the regime, phase, properties, etc. of each region type
        
//         //Dyamical equations: centrifugal and coriolis forces (since using corotating frame)
//     }
// };



// //Gaseous upper atmosphere, with ordinary chemicals and conditions
// struct upper_atmosphere_LAYER : public region {
//     array<double,5> composition_percentages = {10,20.3,10,50,8.7}; //placeholder

//     //EOS and dynamical equations are used for change(). Check if we need to account for phase transitions where, e.g., hydrogen becomes liquid or solid
//     //EOS: ideal gas law with cloud chemistry corrections
//     //DYNAMICS:
//     // compressible Navier-Stokes equations
//     // gravity
//     // radiative heating/cooling
//     // cloud microphysics and condensation/evaporation

//     //add a constructor to each region type w/ all variables to take in (some may be calculated using EOS and therefore not inputted as parameters)
// };

// struct molecular_hydrogen_LAYER : public region {
//     array<double,5> composition_percentages = {1,2,3,4,5}; //placeholder

//     //EOS: non-ideal H-He molecular-fluid EOS (e.g. Saumon-Chabrier-Van Horn / modern ab-initio tables) *Convection should be a natural result of EOS*
//     //DYNAMICS:
//     // compressible Navier-Stokes equations
//     // gravity
//     // radiative heating/cooling
//     // turbulent mixing
// };

// //Supercritical fluid hydrogen, liquid helium
// struct supercritical_fluid_LAYER : public region {
//     //EOS: dense supercritical H-He fluid EOS with pressure ionization and non-ideal interactions
//     //DYNAMICS:
//     // compressible Navier-Stokes equations
//     // gravity
//     // radiative and conductive heat transport
//     // turbulent mixing
//     // phase-separation source terms (helium rain)
// };

// //Liquid metallic hydrogen, plasma, helium rain
// struct core_LAYER : public region {
//     int some_phonon_variable = 67;

//     void interact() override {
//         //for thermal distribution, use blackbody except it is cut off for frequence below plasma frequency threshold
//         //electron plasma impacts path lengths in a certain way
//     }

//     //EOS: liquid metallic-hydrogen EOS with partial electron degeneracy, pressure ionization, composition gradients, heavy-element enrichment, and helium-rain corrections
//     //DYNAMICS:
//     // magnetohydrodynamics (mass, momentum, energy, induction equations)
//     // gravity
//     // conductive heat transport
//     // compositional transport/advection
//     // helium-rain source terms
// };

// //liquid metallic hydrogen, heavy elements, rock, ice
// struct inner_core_LAYER : public region {
    
//     //EOS: metallic-hydrogen + heavy-element mixture EOS with strong electron degeneracy, superionic/dissolved-material corrections, and high-pressure ab-initio equation-of-state tables
//     //DYNAMICS:
//     // magnetohydrodynamics (mass, momentum, energy, induction equations)
//     // gravity
//     // conductive heat transport
//     // compositional transport/advection
//     // heavy-element settling and mixing
// };