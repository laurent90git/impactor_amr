// Command-line entry point: parses the options, builds the model of the chosen case and runs it.
//
//   ./hvi1d --case shear-slab --level 6
//   ./hvi2d --case impact2009 --level 5 --coarse-level 2 --paper-times
//   ./hvi2d --case impact2012 --level 5 --coarse-level 2 --paper-times
//
// The mesh dimension is a compile-time constant (HVI_DIM, one executable per dimension).
#include <iostream>
#include <string>

#include <samurai/samurai.hpp>

#include "../cases/cases.hpp"
#include "../physics/model.hpp"
#include "../solver/run.hpp"

#ifndef HVI_DIM
#define HVI_DIM 2
#endif

int main(int argc, char* argv[])
{
    constexpr std::size_t dim = HVI_DIM;
    auto& app                 = samurai::initialize("Eulerian hyperelastic(-plastic) diffuse-interface solver", argc, argv);

    std::string case_name = dim == 1 ? "shear-slab" : "impact2009";
    double mu             = -1.;
    hvi::RunOptions opt;
    app.add_option("--case", case_name,
                   "advect, slip-advect, shear-slab, double-expansion, solid-shear-tube, solid-gas-tube, impact1d, wilkins, "
                   "impact2009, impact2012, penetration2015")
        ->capture_default_str();
    app.add_option("--level", opt.level, "finest level (default: the case's)");
    app.add_option("--coarse-level", opt.min_level, "coarsest level (< level: multiresolution adaptation)");
    app.add_option("--tf", opt.tf, "final time (default: the case's)");
    app.add_option("--nout", opt.nout, "number of evenly spaced outputs")->capture_default_str();
    app.add_flag("--paper-times", opt.paper_times, "outputs at the reference times of the paper");
    app.add_option("--cfl", opt.cfl, "CFL number (on the sum of the speeds of all directions)")->capture_default_str();
    app.add_option("--limiter", opt.limiter, "none (first order), minmod, vanleer, mc")->capture_default_str();
    app.add_option("--muscl-cobasis", opt.muscl_A, "reconstruct the cobasis (1) or not (0)")->capture_default_str();
    app.add_option("--rk", opt.rk, "1: forward Euler, 2: SSP-RK2")->capture_default_str();
    app.add_option("--alim", opt.alim, "fluidisation threshold on the solid fraction (0: no slip treatment)");
    app.add_option("--plastic", opt.plastic, "1: plasticity on, 0: off (default: the case's)");
    app.add_option("--eps", opt.eps, "residual volume fraction (default: the case's)");
    app.add_option("--mu", mu, "shear modulus of the solid of the 2-phase cases (impact2009: 9.2e10, 1e9, 0)");
    app.add_option("--adapt-eps", opt.mr_eps, "multiresolution threshold")->capture_default_str();
    app.add_option("--adapt-reg", opt.mr_reg, "multiresolution regularity")->capture_default_str();
    app.add_option("--axis", opt.axis, "direction of the 1D cases")->capture_default_str();
    app.add_option("--out", opt.out_dir, "output directory")->capture_default_str();
    app.add_option("--tag", opt.tag, "suffix of the output files");
    app.add_option("--log-every", opt.log_every, "log period in steps")->capture_default_str();
    app.add_option("--max-steps", opt.max_steps, "maximum number of steps");
    SAMURAI_PARSE(argc, argv);

    if (opt.axis >= dim)
    {
        std::cerr << "--axis must be < dim" << std::endl;
        return 1;
    }
    int ret = 0;
    if (hvi::is_3phase_case(case_name))
    {
        auto cd = hvi::make_case_3phase(case_name, opt.axis);
        hvi::Model<3, 2> model(cd.mat);
        ret = hvi::run_case<dim>(model, cd, opt);
    }
    else
    {
        auto cd = hvi::make_case_2phase(case_name, opt.axis, mu);
        hvi::Model<2, 1> model(cd.mat);
        ret = hvi::run_case<dim>(model, cd, opt);
    }
    samurai::finalize();
    return ret;
}
