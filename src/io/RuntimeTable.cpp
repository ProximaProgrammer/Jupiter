#include "RuntimeTable.hpp"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <stdexcept>

namespace Jupiter
{

namespace
{
constexpr std::uint32_t kMagic = 0x4A525431u; // JRT1
constexpr std::uint32_t kVersion = 1u;

template <class T>
void read_exact(std::ifstream& in, T& x)
{
    in.read(reinterpret_cast<char*>(&x), sizeof(T));
    if (!in) throw std::runtime_error("Unexpected end of runtime table");
}

template <class T>
void read_vec(std::ifstream& in, std::vector<T>& v, std::size_t n)
{
    v.resize(n);
    if (n == 0) return;
    in.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(sizeof(T) * n));
    if (!in) throw std::runtime_error("Unexpected end of runtime table vector");
}

std::string read_string(std::ifstream& in)
{
    std::uint32_t n = 0;
    read_exact(in, n);
    std::string s(n, '\0');
    if (n > 0) in.read(s.data(), n);
    if (!in) throw std::runtime_error("Unexpected end of string");
    return s;
}

void check_header(std::ifstream& in)
{
    std::uint32_t magic = 0, version = 0;
    read_exact(in, magic);
    read_exact(in, version);
    if (magic != kMagic || version != kVersion) throw std::runtime_error("Bad runtime-table header/version");
}
} // namespace

bool RuntimeTables::Load(const std::string& asset_dir)
{
    last_error_.clear();
    wavelengths_m.clear();
    opacities_.clear();
    cia_.clear();
    h3plus_ = EmissionTable{};

    try
    {
        {
            std::ifstream in(asset_dir + "/spectral_grid.bin", std::ios::binary);
            if (!in) throw std::runtime_error("Missing spectral_grid.bin in " + asset_dir);
            check_header(in);
            std::uint64_t n = 0;
            read_exact(in, n);
            read_vec(in, wavelengths_m, static_cast<std::size_t>(n));
        }

        {
            std::ifstream in(asset_dir + "/opacities.bin", std::ios::binary);
            if (in)
            {
                check_header(in);
                std::uint32_t nSpecies = 0;
                read_exact(in, nSpecies);
                for (std::uint32_t s = 0; s < nSpecies; ++s)
                {
                    SpeciesOpacityTable tab;
                    tab.name = read_string(in);
                    std::uint64_t nP = 0, nT = 0, nL = 0;
                    read_exact(in, nP); read_exact(in, nT); read_exact(in, nL);
                    read_vec(in, tab.pressures_bar, static_cast<std::size_t>(nP));
                    read_vec(in, tab.temperatures_K, static_cast<std::size_t>(nT));
                    read_vec(in, tab.sigma_m2_per_molecule, static_cast<std::size_t>(nP*nT*nL));
                    opacities_.push_back(std::move(tab));
                }
            }
        }

        {
            std::ifstream in(asset_dir + "/cia.bin", std::ios::binary);
            if (in)
            {
                check_header(in);
                std::uint32_t nPairs = 0;
                read_exact(in, nPairs);
                for (std::uint32_t s = 0; s < nPairs; ++s)
                {
                    CIATable tab;
                    tab.pair_name = read_string(in);
                    std::uint64_t nT = 0, nL = 0;
                    read_exact(in, nT); read_exact(in, nL);
                    read_vec(in, tab.temperatures_K, static_cast<std::size_t>(nT));
                    read_vec(in, tab.coeff_m5_per_molecule2, static_cast<std::size_t>(nT*nL));
                    cia_.push_back(std::move(tab));
                }
            }
        }

        {
            std::ifstream in(asset_dir + "/emission.bin", std::ios::binary);
            if (in)
            {
                check_header(in);
                h3plus_.name = read_string(in);
                std::uint64_t nT = 0, nL = 0;
                read_exact(in, nT); read_exact(in, nL);
                read_vec(in, h3plus_.temperatures_K, static_cast<std::size_t>(nT));
                read_vec(in, h3plus_.F_lambda, static_cast<std::size_t>(nT*nL));
            }
        }
        return true;
    }
    catch (const std::exception& e)
    {
        last_error_ = e.what();
        wavelengths_m.clear();
        opacities_.clear();
        cia_.clear();
        return false;
    }
}

std::size_t RuntimeTables::LowerIndex(const std::vector<double>& axis, double x) const
{
    if (axis.size() < 2) return 0;
    if (x <= axis.front()) return 0;
    if (x >= axis.back()) return axis.size() - 2;
    return static_cast<std::size_t>(std::upper_bound(axis.begin(), axis.end(), x) - axis.begin() - 1);
}

double RuntimeTables::InterpLambda(const std::vector<float>& values, std::size_t offset, std::size_t nLambda, double wavelength_m) const
{
    if (wavelengths_m.size() < 2 || nLambda < 2) return 0.0;
    const std::size_t il = LowerIndex(wavelengths_m, wavelength_m);
    const double x0 = wavelengths_m[il];
    const double x1 = wavelengths_m[il + 1];
    const double f = (x1 > x0) ? (wavelength_m - x0) / (x1 - x0) : 0.0;
    return static_cast<double>(values[offset + il]) * (1.0 - f) + static_cast<double>(values[offset + il + 1]) * f;
}

double RuntimeTables::InterpolateOpacity(const std::string& species, double pressure_bar, double temperature_K, double wavelength_m) const
{
    const auto it = std::find_if(opacities_.begin(), opacities_.end(), [&](const auto& t){ return t.name == species; });
    if (it == opacities_.end() || wavelengths_m.empty()) return 0.0;
    const SpeciesOpacityTable& tab = *it;
    if (tab.pressures_bar.empty() || tab.temperatures_K.empty()) return 0.0;
    const std::size_t nT = tab.temperatures_K.size();
    const std::size_t nL = wavelengths_m.size();
    const std::size_t ip = LowerIndex(tab.pressures_bar, pressure_bar);
    const std::size_t it0 = LowerIndex(tab.temperatures_K, temperature_K);
    const double p0 = tab.pressures_bar[ip], p1 = tab.pressures_bar[ip + 1];
    const double t0 = tab.temperatures_K[it0], t1 = tab.temperatures_K[it0 + 1];
    const double fp = (p1 > p0) ? std::clamp((pressure_bar - p0)/(p1 - p0), 0.0, 1.0) : 0.0;
    const double ft = (t1 > t0) ? std::clamp((temperature_K - t0)/(t1 - t0), 0.0, 1.0) : 0.0;

    const auto sample = [&](std::size_t pidx, std::size_t tidx)
    {
        return InterpLambda(tab.sigma_m2_per_molecule, (pidx*nT + tidx)*nL, nL, wavelength_m);
    };
    const double a00 = sample(ip, it0);
    const double a10 = sample(ip + 1, it0);
    const double a01 = sample(ip, it0 + 1);
    const double a11 = sample(ip + 1, it0 + 1);
    return (1-fp)*(1-ft)*a00 + fp*(1-ft)*a10 + (1-fp)*ft*a01 + fp*ft*a11;
}

double RuntimeTables::InterpolateCIA(const std::string& pair_name, double temperature_K, double wavelength_m) const
{
    const auto it = std::find_if(cia_.begin(), cia_.end(), [&](const auto& t){ return t.pair_name == pair_name; });
    if (it == cia_.end() || wavelengths_m.empty()) return 0.0;
    const CIATable& tab = *it;
    if (tab.temperatures_K.size() < 2) return 0.0;
    const std::size_t nL = wavelengths_m.size();
    const std::size_t it0 = LowerIndex(tab.temperatures_K, temperature_K);
    const double t0 = tab.temperatures_K[it0], t1 = tab.temperatures_K[it0 + 1];
    const double ft = (t1 > t0) ? std::clamp((temperature_K - t0)/(t1 - t0), 0.0, 1.0) : 0.0;
    const double a0 = InterpLambda(tab.coeff_m5_per_molecule2, it0*nL, nL, wavelength_m);
    const double a1 = InterpLambda(tab.coeff_m5_per_molecule2, (it0+1)*nL, nL, wavelength_m);
    return (1.0 - ft)*a0 + ft*a1;
}

double RuntimeTables::InterpolateH3Emission(double temperature_K, double wavelength_m) const
{
    if (h3plus_.temperatures_K.size() < 2 || wavelengths_m.empty()) return 0.0;
    const std::size_t nL = wavelengths_m.size();
    const std::size_t it0 = LowerIndex(h3plus_.temperatures_K, temperature_K);
    const double t0 = h3plus_.temperatures_K[it0], t1 = h3plus_.temperatures_K[it0 + 1];
    const double ft = (t1 > t0) ? std::clamp((temperature_K - t0)/(t1 - t0), 0.0, 1.0) : 0.0;
    const double a0 = InterpLambda(h3plus_.F_lambda, it0*nL, nL, wavelength_m);
    const double a1 = InterpLambda(h3plus_.F_lambda, (it0+1)*nL, nL, wavelength_m);
    return (1.0 - ft)*a0 + ft*a1;
}

} // namespace Jupiter
