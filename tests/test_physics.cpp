// Checks of the physics kernel (no samurai). Run: ./build/test_physics
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "numerics/riemann.hpp"
#include "physics/model.hpp"
#include "physics/relaxation.hpp"

using namespace hvi;

static int failures = 0;

static void check(bool ok, const char* what, double a = 0., double b = 0.)
{
    std::printf("%s %-60s %.10g %.10g\n", ok ? "ok  " : "FAIL", what, a, b);
    if (!ok)
    {
        ++failures;
    }
}

static bool close(double a, double b, double rtol)
{
    return std::abs(a - b) <= rtol * std::max(std::abs(a), std::abs(b)) + 1e-300;
}

int main()
{
    using M = Model<2, 1>;
    const Material cu = materials::copper_2009();
    M model({cu, materials::air()});

    // ---- prim <-> cons round trip ---------------------------------------------------------------------------------
    M::prim_t w;
    w.alpha = {0.7, 0.3};
    w.rho   = {9000., 1.2};
    w.p     = {2e8, 2e8};
    w.v     = {100., -20., 5.};
    w.A[0]  = {{{1.01, 0.02, 0.}, {-0.03, 0.99, 0.01}, {0., 0., 1.}}};
    w.chi   = {1.};
    const auto u  = model.to_cons(w);
    const auto w2 = model.to_prim(u);
    check(close(w2.p[0], w.p[0], 1e-10) && close(w2.p[1], w.p[1], 1e-10), "round trip pressures", w2.p[0], w.p[0]);
    check(close(w2.v[1], w.v[1], 1e-12), "round trip velocity", w2.v[1], w.v[1]);

    // ---- elasticity at rest ---------------------------------------------------------------------------------------
    {
        const auto es = NeoHookeanSeparable::evaluate(1.1 * identity3());
        const Mat3 S  = NeoHookeanSeparable::deviator(es, cu, 8900.);
        check(std::abs(S[0][0]) < 1e-3 && std::abs(NeoHookeanSeparable::energy(es, cu)) < 1e-12, "no deviator for a dilatation", S[0][0]);
    }
    // ---- simple shear: A = I - s e_y (x) e_x  ->  S_01 = mu s (1 + s^2/2) at rho = rho0 --------------------------
    {
        const double s = 0.1;
        Mat3 A         = identity3();
        A[1][0]        = -s;
        const auto es  = NeoHookeanSeparable::evaluate(A);
        const Mat3 S   = NeoHookeanSeparable::deviator(es, cu, cu.rho0);
        check(close(S[0][1], cu.mu * s * (1. + 0.5 * s * s), 1e-12), "simple shear stress", S[0][1], cu.mu * s * (1. + 0.5 * s * s));
    }
    // ---- acoustic tensor at rest: 4/3 mu/rho0 along d, mu/rho0 across ---------------------------------------------
    {
        const auto es = NeoHookeanSeparable::evaluate(identity3());
        const Mat3 D  = NeoHookeanSeparable::acoustic(es, cu, 0);
        check(close(D[0][0], 4. / 3. * cu.mu / cu.rho0, 1e-12) && close(D[1][1], cu.mu / cu.rho0, 1e-12) && std::abs(D[0][1]) < 1e-9,
              "acoustic tensor at rest", D[0][0], 4. / 3. * cu.mu / cu.rho0);
    }
    // ---- uniaxial: D_dd = (mu/rho0)(11/9 a^(8/3) + 1/9 a^(-4/3)) for A = diag(a,1,1) ---------------------------------
    {
        const double a = 1.15;
        Mat3 A         = identity3();
        A[0][0]        = a;
        const auto es  = NeoHookeanSeparable::evaluate(A);
        const Mat3 D   = NeoHookeanSeparable::acoustic(es, cu, 0);
        const double ex = cu.mu / cu.rho0 * (11. / 9. * std::pow(a, 8. / 3.) + 1. / 9. * std::pow(a, -4. / 3.));
        check(close(D[0][0], ex, 1e-12), "uniaxial acoustic tensor", D[0][0], ex);
    }
    // ---- face solver: equal states give the exact flux -------------------------------------------------------------
    {
        LongitudinalHllcShear<M> solver;
        for (std::size_t d = 0; d < 3; ++d)
        {
            const auto r  = solver(model, w, w, d);
            const auto mx = model.mixture(w, d);
            using L       = M::layout_t;
            const double Fm = mx.rho * w.v[d] * w.v[d] - mx.sigma[d][d];
            const double Ft = mx.rho * w.v[d] * w.v[(d + 1) % 3] - mx.sigma[(d + 1) % 3][d];
            const double FE = (mx.rhoE) * w.v[d] - (mx.sigma[0][d] * w.v[0] + mx.sigma[1][d] * w.v[1] + mx.sigma[2][d] * w.v[2]);
            check(close(r.flux[L::mom(d)], Fm, 1e-10) && close(r.flux[L::mom((d + 1) % 3)], Ft, 1e-10) && close(r.flux[L::energy], FE, 1e-10)
                      && close(r.vstar[d], w.v[d], 1e-12),
                  "consistency of the face flux (normal, tangential, energy)", r.flux[L::energy], FE);
        }
    }
    // ---- relaxation leaves an equilibrium state unchanged ------------------------------------------------------------
    {
        CellRelaxation<M> relax;
        auto c = u;
        RelaxationStats st;
        relax(model, c, st);
        const auto w3 = model.to_prim(c);
        check(close(w3.p[0], w.p[0], 1e-9) && close(w3.alpha[0], w.alpha[0], 1e-9), "relaxation of an equilibrium state", w3.p[0], w.p[0]);
    }
    // ---- relaxation of a pressure jump conserves energy, equalises the pressures --------------------------------------
    {
        auto w4 = w;
        w4.p    = {5e8, 1e5};
        auto c  = model.to_cons(w4);
        const double E0 = c[M::layout_t::energy];
        CellRelaxation<M> relax;
        RelaxationStats st;
        relax(model, c, st);
        const auto w5 = model.to_prim(c);
        check(close(w5.p[0], w5.p[1], 1e-9) && close(c[M::layout_t::energy], E0, 1e-12), "pressure relaxation", w5.p[0], w5.p[1]);
    }
    // ---- plasticity: a strongly sheared state ends on the yield surface, density unchanged ---------------------------
    {
        const Material pl = materials::copper_2012();
        Mat3 A            = identity3();
        A[1][0]           = -0.2;
        A[0][0]           = 1.05;
        const double J0   = det(A);
        VonMisesInstantaneous<NeoHookeanSeparable> vm;
        const bool y   = vm(A, pl, 9000.);
        const Mat3 S   = NeoHookeanSeparable::deviator(NeoHookeanSeparable::evaluate(A), pl, 9000.);
        const double q = std::sqrt(1.5 * ddot(S, S));
        check(y && close(q, pl.yield, 1e-6) && close(det(A), J0, 1e-12), "von Mises return (equivalent stress = yield)", q, pl.yield);
    }

    std::printf("%d failure(s)\n", failures);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
