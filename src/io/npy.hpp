// Minimal .npy writer (float64, C order).
#pragma once

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace hvi
{
    inline void write_npy(const std::string& path, const std::vector<double>& a, const std::vector<std::size_t>& shape)
    {
        std::string dict = "{'descr': '<f8', 'fortran_order': False, 'shape': (";
        for (auto n : shape)
        {
            dict += std::to_string(n) + ", ";
        }
        dict += "), }";
        while ((10 + dict.size() + 1) % 64 != 0)
        {
            dict += ' ';
        }
        dict += '\n';
        std::ofstream f(path, std::ios::binary);
        const char magic[] = "\x93NUMPY\x01\x00";
        f.write(magic, 8);
        const auto len = static_cast<std::uint16_t>(dict.size());
        f.write(reinterpret_cast<const char*>(&len), 2);
        f << dict;
        f.write(reinterpret_cast<const char*>(a.data()), static_cast<std::streamsize>(a.size() * sizeof(double)));
    }
} // namespace hvi
