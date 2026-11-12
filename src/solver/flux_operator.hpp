// samurai glue: the space operator of the hyperbolic step, as a samurai flux-based scheme.
//
// For each face j+1/2 in direction d (stencil of 4 cells j-1 .. j+2):
//   primitive states of the 4 cells -> reconstruction (left / right face states) -> face solver (flux F, velocity v*).
// The scheme returns, in each cell, the divergence  (contribution of the right face + contribution of the left face)
// / h  (the caller writes  U - dt * op(U)). With samurai's non-conservative flux pairs:
//   contribution to the LEFT cell of the face (flux[0])  =   F + NC(cell j,   v*)
//   contribution to the RIGHT cell of the face (flux[1]) = -(F + NC(cell j+1, v*))
// so that a cell i receives  F_{i+1/2} - F_{i-1/2} + NC(i, v*_{i+1/2} - v*_{i-1/2}), NC being linear in v*.
// NC holds the non-conservative products of the model, evaluated with the cell's own values:
//   d alpha_l / dt          += alpha_l       d u*_d / dx_d
//   d (alpha rho e)_l / dt  -= alpha_l p_l   d u*_d / dx_d
//   d A[b][i] / dt          += A[b][i]       d u*_d / dx_d        (i != d)
//   d A[b][d] / dt          -= sum_{t != d} A[b][t] d v*_t / dx_d
// (multi-dimensional form of  D e^b / Dt + (grad v)^T e^b = 0, Ndanou 2015 eq. 15, summed over the directions).
#pragma once

#include <array>
#include <cstddef>

#include <samurai/schemes/fv.hpp>

#include "../numerics/reconstruction.hpp"
#include "../numerics/riemann.hpp"

namespace hvi
{
    /// NC(cell, v): see the header comment. Adds the non-conservative contribution of the face velocity v to flux.
    template <class Model>
    void add_nonconservative(const Model& model, const typename Model::prim_t& w, const Vec3& v, std::size_t d, typename Model::cons_t& f)
    {
        using L = typename Model::layout_t;
        for (std::size_t l = 0; l < Model::N; ++l)
        {
            f[L::alpha(l)] -= w.alpha[l] * v[d];
            f[L::arhoe(l)] += w.alpha[l] * w.p[l] * v[d];
        }
        (void)model;
        for (std::size_t s = 0; s < Model::NS; ++s)
        {
            for (std::size_t b = 0; b < 3; ++b)
            {
                for (std::size_t t = 0; t < 3; ++t)
                {
                    if (t == d)
                    {
                        continue;
                    }
                    f[L::A(s, b, t)] -= w.A[s][b][t] * v[d];
                    f[L::A(s, b, d)] += w.A[s][b][t] * v[t];
                }
            }
        }
    }

    /// Builds the samurai operator. Model, Reconstruction and FaceSolver are copied into the lambdas.
    template <class Field, class Model, class Reconstruction = MusclPrimitive<Model>, class FaceSolver = LongitudinalHllcShear<Model>>
    auto make_flux_operator(const Model& model, const Reconstruction& recon = {}, const FaceSolver& solver = {})
    {
        static constexpr std::size_t dim = Field::dim;
        static constexpr std::size_t nc  = Model::nc;
        static_assert(Field::n_comp == nc, "the field must have the components of the model layout");
        using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, 4, Field, Field>;

        samurai::FluxDefinition<cfg> def;
        samurai::static_for<0, dim>::apply(
            [&](auto integral_constant_d)
            {
                static constexpr std::size_t d = decltype(integral_constant_d)::value;
                def[d].flux_function = [model, recon, solver](samurai::FluxValuePair<cfg>& flux,
                                                              const samurai::StencilData<cfg>&,
                                                              const samurai::StencilValues<cfg>& u)
                {
                    using prim_t = typename Model::prim_t;
                    std::array<prim_t, 4> w;
                    for (std::size_t k = 0; k < 4; ++k)
                    {
                        std::array<double, nc> c;
                        for (std::size_t q = 0; q < nc; ++q)
                        {
                            c[q] = u[k](q);
                        }
                        w[k] = model.to_prim(c);
                    }
                    prim_t wl, wr;
                    recon(model, w, wl, wr);
                    const auto res = solver(model, wl, wr, d);
                    auto fl        = res.flux;
                    auto fr        = res.flux;
                    add_nonconservative(model, w[1], res.vstar, d, fl);
                    add_nonconservative(model, w[2], res.vstar, d, fr);
                    for (std::size_t q = 0; q < nc; ++q)
                    {
                        flux[0](q) = fl[q];
                        flux[1](q) = -fr[q];
                    }
                };
            });
        auto scheme = samurai::make_flux_based_scheme(def);
        scheme.set_name("hvi_hyperbolic");
        return scheme;
    }
} // namespace hvi
