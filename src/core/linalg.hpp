// Small fixed-size linear algebra (3-vectors, 3x3 matrices) used by the physics kernel.
// No dependency on samurai: the kernel can be tested and reused on its own.
#pragma once

#include <array>
#include <cmath>
#include <utility>

namespace hvi
{
    using Vec3 = std::array<double, 3>;
    using Mat3 = std::array<std::array<double, 3>, 3>; // m[row][col]

    inline constexpr Mat3 identity3()
    {
        return {{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}};
    }

    inline Mat3 zero3()
    {
        return {};
    }

    inline Mat3 operator*(const Mat3& a, const Mat3& b)
    {
        Mat3 c{};
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
            }
        }
        return c;
    }

    inline Mat3 operator*(double s, const Mat3& a)
    {
        Mat3 c = a;
        for (auto& r : c)
        {
            for (auto& x : r)
            {
                x *= s;
            }
        }
        return c;
    }

    inline Mat3 operator+(const Mat3& a, const Mat3& b)
    {
        Mat3 c = a;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                c[i][j] += b[i][j];
            }
        }
        return c;
    }

    inline Mat3& operator+=(Mat3& a, const Mat3& b)
    {
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                a[i][j] += b[i][j];
            }
        }
        return a;
    }

    inline Mat3 transpose(const Mat3& a)
    {
        Mat3 t{};
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                t[i][j] = a[j][i];
            }
        }
        return t;
    }

    /// a^T a
    inline Mat3 ata(const Mat3& a)
    {
        Mat3 c{};
        for (int i = 0; i < 3; ++i)
        {
            for (int j = i; j < 3; ++j)
            {
                c[i][j] = a[0][i] * a[0][j] + a[1][i] * a[1][j] + a[2][i] * a[2][j];
                c[j][i] = c[i][j];
            }
        }
        return c;
    }

    inline double det(const Mat3& a)
    {
        return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
             + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    }

    inline double trace(const Mat3& a)
    {
        return a[0][0] + a[1][1] + a[2][2];
    }

    /// a : b
    inline double ddot(const Mat3& a, const Mat3& b)
    {
        double s = 0.;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                s += a[i][j] * b[i][j];
            }
        }
        return s;
    }

    /// Eigen-decomposition of a symmetric matrix by cyclic Jacobi rotations: a = V diag(lambda) V^T, the eigenvectors
    /// are the COLUMNS of V.
    inline void eigen_sym(const Mat3& a_in, Vec3& lambda, Mat3& V)
    {
        Mat3 a = a_in;
        V      = identity3();
        for (int sweep = 0; sweep < 50; ++sweep)
        {
            const double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
            const double dia = a[0][0] * a[0][0] + a[1][1] * a[1][1] + a[2][2] * a[2][2];
            if (off <= 1e-30 * dia || off == 0.)
            {
                break;
            }
            for (int p = 0; p < 2; ++p)
            {
                for (int q = p + 1; q < 3; ++q)
                {
                    if (a[p][q] == 0.)
                    {
                        continue;
                    }
                    const double theta = (a[q][q] - a[p][p]) / (2. * a[p][q]);
                    const double t     = (theta >= 0. ? 1. : -1.) / (std::abs(theta) + std::sqrt(theta * theta + 1.));
                    const double c     = 1. / std::sqrt(t * t + 1.);
                    const double s     = t * c;
                    for (int k = 0; k < 3; ++k) // a <- J^T a J
                    {
                        const double akp = a[k][p], akq = a[k][q];
                        a[k][p]          = c * akp - s * akq;
                        a[k][q]          = s * akp + c * akq;
                    }
                    for (int k = 0; k < 3; ++k)
                    {
                        const double apk = a[p][k], aqk = a[q][k];
                        a[p][k]          = c * apk - s * aqk;
                        a[q][k]          = s * apk + c * aqk;
                    }
                    for (int k = 0; k < 3; ++k)
                    {
                        const double vkp = V[k][p], vkq = V[k][q];
                        V[k][p]          = c * vkp - s * vkq;
                        V[k][q]          = s * vkp + c * vkq;
                    }
                }
            }
        }
        lambda = {a[0][0], a[1][1], a[2][2]};
    }

    inline double sq(double x)
    {
        return x * x;
    }
} // namespace hvi
