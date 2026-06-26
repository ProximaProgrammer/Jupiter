#include "io/Visualizer.h"
#include "physics/RadiativeTransfer.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace jrt {
namespace {

struct RGB { unsigned char r=0, g=0, b=0; };

class Image {
public:
    int w, h;
    std::vector<RGB> pix;
    Image(int width, int height) : w(width), h(height), pix(static_cast<std::size_t>(w*h), RGB{255,255,255}) {}
    void set(int x, int y, RGB c) {
        if (x < 0 || y < 0 || x >= w || y >= h) return;
        pix[static_cast<std::size_t>(y*w + x)] = c;
    }
    void rect(int x0, int y0, int x1, int y1, RGB c) {
        for (int y = y0; y < y1; ++y) for (int x = x0; x < x1; ++x) set(x,y,c);
    }
    void writePpm(const std::string& path) const {
        std::ofstream f(path, std::ios::binary);
        if (!f) throw std::runtime_error("Could not write " + path);
        f << "P6\n" << w << ' ' << h << "\n255\n";
        for (const auto& p : pix) { f.put(static_cast<char>(p.r)); f.put(static_cast<char>(p.g)); f.put(static_cast<char>(p.b)); }
    }
};

const std::array<unsigned char, 7>& glyph(char ch) {
    static const std::map<char, std::array<unsigned char, 7>> font = {
        {'A',{0x0E,0x11,0x11,0x1F,0x11,0x11,0x11}}, {'B',{0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E}},
        {'C',{0x0E,0x11,0x10,0x10,0x10,0x11,0x0E}}, {'D',{0x1E,0x11,0x11,0x11,0x11,0x11,0x1E}},
        {'E',{0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F}}, {'F',{0x1F,0x10,0x10,0x1E,0x10,0x10,0x10}},
        {'G',{0x0E,0x11,0x10,0x17,0x11,0x11,0x0E}}, {'H',{0x11,0x11,0x11,0x1F,0x11,0x11,0x11}},
        {'I',{0x0E,0x04,0x04,0x04,0x04,0x04,0x0E}}, {'J',{0x07,0x02,0x02,0x02,0x12,0x12,0x0C}},
        {'K',{0x11,0x12,0x14,0x18,0x14,0x12,0x11}}, {'L',{0x10,0x10,0x10,0x10,0x10,0x10,0x1F}},
        {'M',{0x11,0x1B,0x15,0x15,0x11,0x11,0x11}}, {'N',{0x11,0x19,0x15,0x13,0x11,0x11,0x11}},
        {'O',{0x0E,0x11,0x11,0x11,0x11,0x11,0x0E}}, {'P',{0x1E,0x11,0x11,0x1E,0x10,0x10,0x10}},
        {'Q',{0x0E,0x11,0x11,0x11,0x15,0x12,0x0D}}, {'R',{0x1E,0x11,0x11,0x1E,0x14,0x12,0x11}},
        {'S',{0x0F,0x10,0x10,0x0E,0x01,0x01,0x1E}}, {'T',{0x1F,0x04,0x04,0x04,0x04,0x04,0x04}},
        {'U',{0x11,0x11,0x11,0x11,0x11,0x11,0x0E}}, {'V',{0x11,0x11,0x11,0x11,0x0A,0x0A,0x04}},
        {'W',{0x11,0x11,0x11,0x15,0x15,0x1B,0x11}}, {'X',{0x11,0x0A,0x04,0x04,0x04,0x0A,0x11}},
        {'Y',{0x11,0x0A,0x04,0x04,0x04,0x04,0x04}}, {'Z',{0x1F,0x01,0x02,0x04,0x08,0x10,0x1F}},
        {'0',{0x0E,0x11,0x13,0x15,0x19,0x11,0x0E}}, {'1',{0x04,0x0C,0x04,0x04,0x04,0x04,0x0E}},
        {'2',{0x0E,0x11,0x01,0x02,0x04,0x08,0x1F}}, {'3',{0x1E,0x01,0x01,0x0E,0x01,0x01,0x1E}},
        {'4',{0x02,0x06,0x0A,0x12,0x1F,0x02,0x02}}, {'5',{0x1F,0x10,0x10,0x1E,0x01,0x01,0x1E}},
        {'6',{0x0E,0x10,0x10,0x1E,0x11,0x11,0x0E}}, {'7',{0x1F,0x01,0x02,0x04,0x08,0x08,0x08}},
        {'8',{0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E}}, {'9',{0x0E,0x11,0x11,0x0F,0x01,0x01,0x0E}},
        {' ',{0,0,0,0,0,0,0}}, {'-',{0,0,0,0x1F,0,0,0}}, {'_',{0,0,0,0,0,0,0x1F}}, {'.',{0,0,0,0,0,0x0C,0x0C}},
        {':',{0,0x0C,0x0C,0,0x0C,0x0C,0}}, {'/',{0x01,0x02,0x02,0x04,0x08,0x08,0x10}}, {'=',{0,0,0x1F,0,0x1F,0,0}},
        {'[',{0x0E,0x08,0x08,0x08,0x08,0x08,0x0E}}, {']',{0x0E,0x02,0x02,0x02,0x02,0x02,0x0E}}, {'+',{0,0x04,0x04,0x1F,0x04,0x04,0}},
        {'(',{0x02,0x04,0x08,0x08,0x08,0x04,0x02}}, {')',{0x08,0x04,0x02,0x02,0x02,0x04,0x08}}, {'%',{0x18,0x19,0x02,0x04,0x08,0x13,0x03}}
    };
    static const std::array<unsigned char,7> unknown = {0,0,0x1F,0x01,0x0E,0,0x04};
    const char up = (ch >= 'a' && ch <= 'z') ? static_cast<char>(ch - 'a' + 'A') : ch;
    auto it = font.find(up);
    return it == font.end() ? unknown : it->second;
}

void text(Image& im, int x, int y, const std::string& s, RGB c, int scale=2) {
    int cx = x;
    for (char ch : s) {
        const auto& g = glyph(ch);
        for (int row = 0; row < 7; ++row) {
            for (int col = 0; col < 5; ++col) {
                if (g[static_cast<std::size_t>(row)] & (1 << (4-col))) {
                    for (int yy=0; yy<scale; ++yy) for (int xx=0; xx<scale; ++xx) im.set(cx + col*scale + xx, y + row*scale + yy, c);
                }
            }
        }
        cx += 6 * scale;
    }
}

RGB cmap(double x) {
    const double q = std::clamp(x, 0.0, 1.0);
    const auto b = [](double v) { return static_cast<unsigned char>(std::clamp(v, 0.0, 255.0)); };
    return {b(255.0*q), b(255.0*(1.0-std::abs(q-0.5)*1.3)), b(255.0*(1.0-q))};
}

std::string num(double x, int prec=2) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(prec) << x;
    return os.str();
}

template <class F>
void panelMap(Image& im, const Grid& g, int x0, int y0, int pw, int ph, const std::string& title, F valueFn) {
    const std::size_t i = g.outerInteriorI();
    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -vmin;
    for (std::size_t j = g.jBegin(); j < g.jEnd(); ++j) {
        for (std::size_t k = g.kBegin(); k < g.kEnd(); ++k) {
            const double v = valueFn(i,j,k);
            vmin = std::min(vmin, v); vmax = std::max(vmax, v);
        }
    }
    im.rect(x0, y0, x0+pw, y0+ph, RGB{245,245,245});
    for (int y = 0; y < ph-36; ++y) {
        const double jf = static_cast<double>(y) / std::max(1, ph-37);
        std::size_t j = g.jBegin() + std::min<std::size_t>(g.cfg.ntheta-1, static_cast<std::size_t>(jf * g.cfg.ntheta));
        for (int x = 0; x < pw; ++x) {
            const double kf = static_cast<double>(x) / std::max(1, pw-1);
            std::size_t k = std::min<std::size_t>(g.cfg.nphi-1, static_cast<std::size_t>(kf * g.cfg.nphi));
            const double v = valueFn(i,j,k);
            im.set(x0+x, y0+24+y, cmap((v-vmin) / std::max(1.0e-30, vmax-vmin)));
        }
    }
    text(im, x0+6, y0+5, title, RGB{0,0,0}, 2);
    text(im, x0+6, y0+ph-18, "MIN=" + num(vmin,2) + " MAX=" + num(vmax,2), RGB{0,0,0}, 1);
}

} // namespace

void writeDashboard(const Grid& grid, int step, double time_s, const std::string& outputDir) {
    std::filesystem::create_directories(outputDir + "/frames");
    Image im(1200, 820);
    im.rect(0,0,1200,820,RGB{235,235,235});

    text(im, 20, 18, "JUPITERRT PATCH 04: OUTER INTERIOR SHELL (OUTER RADIAL LAYER IS GHOST)", RGB{0,0,0}, 2);
    text(im, 20, 46, "STEP=" + std::to_string(step) + " TIME_S=" + num(time_s,1) + " SHELL_I=" + std::to_string(grid.outerInteriorI()), RGB{0,0,0}, 2);

    const int pw = 560, ph = 345;
    panelMap(im, grid, 20, 90, pw, ph, "T [K]", [&](std::size_t i, std::size_t j, std::size_t k){ return grid.at(i,j,k).T; });
    panelMap(im, grid, 620, 90, pw, ph, "U_PHI COROT [M/S]", [&](std::size_t i, std::size_t j, std::size_t k){ return grid.at(i,j,k).u.ph; });
    panelMap(im, grid, 20, 455, pw, ph, "U_R [M/S]", [&](std::size_t i, std::size_t j, std::size_t k){ return grid.at(i,j,k).u.r; });
    panelMap(im, grid, 620, 455, pw, ph, "KAPPA 4.8UM [M2/KG]", [&](std::size_t i, std::size_t j, std::size_t k){ return spectralOpacityM2PerKg(grid.at(i,j,k), 4.8e-6); });

    std::ostringstream name;
    name << outputDir << "/frames/dashboard_" << std::setw(6) << std::setfill('0') << step << ".ppm";
    im.writePpm(name.str());
    std::filesystem::copy_file(name.str(), outputDir + "/frames/dashboard_latest.ppm", std::filesystem::copy_options::overwrite_existing);

    std::ofstream meta(name.str() + ".txt");
    meta << "JupiterRT Patch 04 dashboard\n"
         << "The rendered shell is radial index " << grid.outerInteriorI() << "; radial index " << grid.outerGhostI() << " is the outer ghost boundary.\n"
         << "Panels: temperature, co-rotating zonal velocity, radial velocity, 4.8 micron opacity proxy.\n"
         << "Latitude runs north-to-south vertically; longitude runs left-to-right.\n";
}

} // namespace jrt
