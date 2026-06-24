#pragma once

namespace Jupiter
{

struct PrimitiveState
{
    double density = 0.0;

    double pressure = 0.0;

    double temperature = 0.0;

    double velocity_r = 0.0;
    double velocity_theta = 0.0;
    double velocity_phi = 0.0;

    double magnetic_r = 0.0;
    double magnetic_theta = 0.0;
    double magnetic_phi = 0.0;
};

}