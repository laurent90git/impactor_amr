// Second-order reconstruction of the primitive variables at a face (MUSCL), from the 4 cells j-1, j, j+1, j+2 of the
// face j+1/2. To bring another reconstruction, write a struct with the same operator() (stencil of primitive states
// -> left and right face states) and pass it to the flux operator (flux_operator.hpp).
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <stdexcept>

namespace hvi
{
    enum class Limiter
    {
        none,   ///< first order (zero slopes)
        minmod,
        vanleer,
        mc      ///< monotonised central
    };

    inline Limiter limiter_from_string(const std::string& s)
    {
        if (s == "none" || s == "o1")
        {
            return Limiter::none;
        }
        if (s == "minmod")
        {
            return Limiter::minmod;
        }
        if (s == "vanleer")
        {
            return Limiter::vanleer;
        }
        if (s == "mc")
        {
            return Limiter::mc;
        }
        throw std::invalid_argument("unknown limiter " + s);
    }

    /// limited slope from the left and right differences
    inline double limited_slope(double dm, double dp, Limiter lim)
    {
        if (dm * dp <= 0.)
        {
            return 0.;
        }
        switch (lim)
        {
            case Limiter::none: return 0.;
            case Limiter::minmod: return std::abs(dm) < std::abs(dp) ? dm : dp;
            case Limiter::vanleer: return 2. * dm * dp / (dm + dp);
            case Limiter::mc:
            {
                const double s = dm > 0. ? 1. : -1.;
                return s * std::min({2. * std::abs(dm), 2. * std::abs(dp), 0.5 * std::abs(dm + dp)});
            }
        }
        return 0.;
    }

    /// MUSCL on the primitive variables (alpha_l, rho_l, v, p_l and optionally the cobasis). chi is not reconstructed:
    /// each face state keeps the chi of its cell. A face state that is not admissible (negative density, pressure
    /// below -p_inf) falls back to the cell value.
    template <class Model>
    struct MusclPrimitive
    {
        using prim_t = typename Model::prim_t;
        Limiter limiter      = Limiter::minmod;
        bool reconstruct_A   = true;

        void operator()(const Model& model, const std::array<prim_t, 4>& w, prim_t& wl, prim_t& wr) const
        {
            wl = w[1];
            wr = w[2];
            if (limiter == Limiter::none)
            {
                return;
            }
            // value at the right face of cell c (sign +1) or at its left face (sign -1)
            auto edge = [&](double a, double b, double c) { return limited_slope(b - a, c - b, limiter); };
            constexpr std::size_t N = Model::N, NS = Model::NS;
            for (std::size_t l = 0; l < N; ++l)
            {
                wl.alpha[l] = w[1].alpha[l] + 0.5 * edge(w[0].alpha[l], w[1].alpha[l], w[2].alpha[l]);
                wr.alpha[l] = w[2].alpha[l] - 0.5 * edge(w[1].alpha[l], w[2].alpha[l], w[3].alpha[l]);
                wl.rho[l]   = w[1].rho[l] + 0.5 * edge(w[0].rho[l], w[1].rho[l], w[2].rho[l]);
                wr.rho[l]   = w[2].rho[l] - 0.5 * edge(w[1].rho[l], w[2].rho[l], w[3].rho[l]);
                wl.p[l]     = w[1].p[l] + 0.5 * edge(w[0].p[l], w[1].p[l], w[2].p[l]);
                wr.p[l]     = w[2].p[l] - 0.5 * edge(w[1].p[l], w[2].p[l], w[3].p[l]);
            }
            for (std::size_t i = 0; i < 3; ++i)
            {
                wl.v[i] = w[1].v[i] + 0.5 * edge(w[0].v[i], w[1].v[i], w[2].v[i]);
                wr.v[i] = w[2].v[i] - 0.5 * edge(w[1].v[i], w[2].v[i], w[3].v[i]);
            }
            if (reconstruct_A)
            {
                for (std::size_t s = 0; s < NS; ++s)
                {
                    for (std::size_t b = 0; b < 3; ++b)
                    {
                        for (std::size_t i = 0; i < 3; ++i)
                        {
                            wl.A[s][b][i] = w[1].A[s][b][i] + 0.5 * edge(w[0].A[s][b][i], w[1].A[s][b][i], w[2].A[s][b][i]);
                            wr.A[s][b][i] = w[2].A[s][b][i] - 0.5 * edge(w[1].A[s][b][i], w[2].A[s][b][i], w[3].A[s][b][i]);
                        }
                    }
                }
            }
            sanitize(model, wl, w[1]);
            sanitize(model, wr, w[2]);
        }

      private:

        static void sanitize(const Model& model, prim_t& f, const prim_t& cell)
        {
            double sa = 0.;
            bool ok   = true;
            for (std::size_t l = 0; l < Model::N; ++l)
            {
                ok = ok && f.alpha[l] > 0. && f.rho[l] > 0. && f.p[l] > model.mat[l].pmin();
                sa += f.alpha[l];
            }
            if (!ok)
            {
                f = cell;
                return;
            }
            for (auto& a : f.alpha)
            {
                a /= sa;
            }
        }
    };
} // namespace hvi
