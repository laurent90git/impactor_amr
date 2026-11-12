// Boundary conditions filling all the ghost layers at once (samurai's Neumann condition fills only the first one).
// Stencil of a condition: cells[0 .. R-1] inside (cells[R-1] next to the boundary), cells[R .. 2R-1] ghosts, outwards.
#pragma once

#include <samurai/bc.hpp>

namespace hvi
{
    /// zero-order extrapolation of the boundary cell (transmissive / outflow)
    template <class Field, std::size_t R>
    struct TransmissiveImpl : public samurai::Bc<Field>
    {
        INIT_BC(TransmissiveImpl, 2 * R)

        stencil_t get_stencil(constant_stencil_size_t) const override
        {
            return samurai::line_stencil<dim, 0, 2 * R>();
        }

        apply_function_t get_apply_function(constant_stencil_size_t, const direction_t&) const override
        {
            return [](Field& f, const stencil_cells_t& cells, const value_t&)
            {
                for (std::size_t k = R; k < 2 * R; ++k)
                {
                    f[cells[k]] = f[cells[R - 1]];
                }
            };
        }
    };

    template <std::size_t R>
    struct TransmissiveR
    {
        template <class Field>
        using impl_t = TransmissiveImpl<Field, R>;
    };

    using Transmissive = TransmissiveR<2>;
} // namespace hvi
