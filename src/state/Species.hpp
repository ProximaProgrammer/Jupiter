#pragma once

#include <array>
#include <cstddef>
#include <string_view>

namespace Jupiter
{

enum class Species : std::size_t
{
    H2 = 0,
    He,
    CH4,
    NH3,
    PH3,
    H2O,
    H3plus,
    Count
};

inline constexpr std::size_t NumSpecies = static_cast<std::size_t>(Species::Count);

inline constexpr std::array<std::string_view, NumSpecies> SpeciesNames = {
    "H2", "He", "CH4", "NH3", "PH3", "H2O", "H3plus"
};

inline constexpr std::array<double, NumSpecies> SpeciesMolarMassKgPerMol = {
    2.01588e-3,
    4.002602e-3,
    16.04246e-3,
    17.03052e-3,
    33.99758e-3,
    18.01528e-3,
    3.02382e-3
};

inline std::size_t SIndex(Species s) { return static_cast<std::size_t>(s); }

} // namespace Jupiter
