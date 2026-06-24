#pragma once

namespace Jupiter
{

struct CellGeometry
{
    double r_inner = 0.0;
    double r_outer = 0.0;

    double theta_lower = 0.0;
    double theta_upper = 0.0;

    double phi_lower = 0.0;
    double phi_upper = 0.0;

    double volume = 0.0;

    double area_r_minus = 0.0;
    double area_r_plus = 0.0;

    double area_theta_minus = 0.0;
    double area_theta_plus = 0.0;

    double area_phi_minus = 0.0;
    double area_phi_plus = 0.0;
};

}