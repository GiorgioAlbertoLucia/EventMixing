#pragma once

#include <array>
#include <cmath>
#include <TRandom3.h>

namespace physics {

    namespace mass {
        const float kMassPion = 0.13957; // GeV/c^2
        const float kProton = 0.938272; // GeV/c^2
        const float kNeutron = 0.939565; // GeV/c^2
        const float kHelium3 = 2.8083916; // GeV/c^2
    }

    /**
     * Calculate the invariant mass of two particles given their momentum vectors and masses.
     * @param p1 Momentum vector of the first particle (p, eta, phi).
     * @param p2 Momentum vector of the second particle (p, eta, phi).
     * @param m1 Mass of the first particle.
     * @param m2 Mass of the second particle.
    */
    float invariantMass(const std::array<float, 3>& p1, const std::array<float, 3>& p2, float m1, float m2) {
        float px1 = p1[0] * std::cos(p1[2]);
        float py1 = p1[0] * std::sin(p1[2]);
        float pz1 = p1[0] * std::sinh(p1[1]);
        float e1 = sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + m1 * m1);
        float px2 = p2[0] * std::cos(p2[2]);
        float py2 = p2[0] * std::sin(p2[2]);
        float pz2 = p2[0] * std::sinh(p2[1]);
        float e2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + m2 * m2);
        float eTot = e1 + e2;
        float pxTot = px1 + px2;
        float pyTot = py1 + py2;
        float pzTot = pz1 + pz2;
        return std::sqrt(eTot * eTot - (pxTot * pxTot + pyTot * pyTot + pzTot * pzTot));
    }

    float momentumMother(const std::array<float, 3>& p1, const std::array<float, 3>& p2) {
        float px1 = p1[0] * cos(p1[2]);
        float py1 = p1[0] * sin(p1[2]);
        float pz1 = p1[0] * sinh(p1[1]);
        float px2 = p2[0] * cos(p2[2]);
        float py2 = p2[0] * sin(p2[2]);
        float pz2 = p2[0] * sinh(p2[1]);
        return sqrt((px1 + px2) * (px1 + px2) + (py1 + py2) * (py1 + py2) + (pz1 + pz2) * (pz1 + pz2));
    }

    float randomAngleRotation(const float phi) {
        float randomAngle = gRandom->Uniform(0, 2 * M_PI);
        if (phi + randomAngle > M_PI) {
            return phi + randomAngle - 2 * M_PI;
        } else if (phi + randomAngle < -M_PI) {
            return phi + randomAngle + 2 * M_PI;
        }
        return phi + randomAngle;
    }

} // namespace physics

