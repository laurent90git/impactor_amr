"""1D reference prototype of the multi-material hyperelastic diffuse-interface model.

Model: Ndanou, Favrie, Gavrilyuk (JCP 295, 2015), elastic part only (no plasticity yet):
  - N phases, one velocity, pressure non-equilibrium + instantaneous pressure relaxation (Saurel et al. 2009);
  - each phase l carries A_l = F_l^{-1} (rows = cobasis vectors e^beta = grad X^beta); stress from g = G/|G|^(1/3),
    G = A^T A, shear energy e^e = mu/(8 rho0) (tr(g^2) - 3)  (Ndanou 2015 eq. 5).
Numerics (differs from the paper on purpose, see DEVELOPMENT_PLAN.md):
  - unsplit sum of a longitudinal HLLC flux (tangential momentum only advected) and a shear two-wave solver using shear
    impedances (after Favrie, Gavrilyuk, Ndanou, JCP 270, 2014, sub-systems 1 and 2/3);
  - "fluidisation": where the total solid volume fraction is below alpha_lim, the solid has no shear strength (its
    cobasis is reset to a pure dilatation, shear impedance 0). This gives slip at solid/fluid interfaces with no level
    set and no ghost fluid;
  - MUSCL (minmod, primitive variables) + SSP-RK2, relaxation after each stage; or first order.

The state is always 3D (3 velocity components, 3x3 A); x is the only space direction here.

Usage:  python hvi1d.py <case> [--n 1000] [--order 2] [--no-slip] [--alim 0.5] [--cfl 0.4] [--plot out.png]
Cases:  advect, slip-advect (Favrie 2009 sec. 6.1), shear-slab (sec. 6.2), solid-shear-tube (sec. 5.5),
        solid-gas-tube (sec. 5.3), double-expansion (sec. 5.4), impact1d
"""

import argparse
import sys

import numpy as np

I3 = np.eye(3)


class Material:
    def __init__(self, name, gamma, pinf, mu=0.0, rho0=1.0):
        self.name, self.gamma, self.pinf, self.mu, self.rho0 = name, gamma, pinf, mu, rho0


COPPER = Material("copper", 4.22, 342.0e8, 9.2e10, 8900.0)
AIR = Material("air", 1.4, 0.0, 0.0, 1.0)


# ----------------------------------------------------------------------------------------------------------------------
# constitutive relations
# ----------------------------------------------------------------------------------------------------------------------
def elastic(A, mu, rho0):
    """A (...,3,3) -> g, g^2, j2 = tr(g^2), specific shear energy e^e, and S/rho_l (deviatoric stress per unit density)."""
    G = np.einsum("...bi,...bj->...ij", A, A)
    J = np.abs(np.linalg.det(A))
    g = G / (J ** (2.0 / 3.0))[..., None, None]
    g2 = g @ g
    j2 = np.trace(g2, axis1=-2, axis2=-1)
    ee = mu / (8.0 * rho0) * (j2 - 3.0)
    S_over_rho = -(mu / (2.0 * rho0)) * (g2 - (j2 / 3.0)[..., None, None] * I3)
    return g, g2, j2, ee, S_over_rho


def acoustic_elastic(g, g2, j2, mu, rho0, d=0):
    """Elastic part D (per unit mass) of the acoustic tensor in direction d, for e^e = mu/(8 rho0)(j2 - 3):
    D = mu/(2 rho0) [ g_dd g + (g e_d)(g e_d)^T + g^2 - ((g^2 e_d) e_d^T + e_d (g^2 e_d)^T)/3 + j2/9 e_d e_d^T ]."""
    ed = I3[d]
    gd = g[..., :, d]
    g2d = g2[..., :, d]
    D = (
        g[..., d, d][..., None, None] * g
        + gd[..., :, None] * gd[..., None, :]
        + g2
        - (g2d[..., :, None] * ed[None, :] + ed[:, None] * g2d[..., None, :]) / 3.0
        + (j2 / 9.0)[..., None, None] * np.outer(ed, ed)
    )
    return (mu / (2.0 * rho0)) * D


class Model:
    def __init__(self, mats, alim=0.5, slip=True, floor=True):
        self.mats = mats
        self.N = len(mats)
        self.gam = np.array([m.gamma for m in mats])[:, None]
        self.pinf = np.array([m.pinf for m in mats])[:, None]
        self.solid = np.array([m.mu > 0.0 for m in mats])
        self.alim = alim
        self.slip = slip  # False: the solid keeps its shear strength whatever its volume fraction (no slip treatment)
        self.floor = floor  # True: relaxed pressure kept where the energy-based pressure is not admissible
        self.nfloor = 0

    # -- state: dict al (N,n), ar (N,n) = alpha rho, m (3,n) = rho v, E (n) = rho E, ae (N,n) = alpha rho e^h, A (N,n,3,3)
    def chi(self, al):
        if not self.slip:
            return np.ones(al.shape[1])
        return (al[self.solid].sum(axis=0) >= self.alim).astype(float)

    def primitives(self, U):
        al = U["al"]
        rl = U["ar"] / al
        rho = U["ar"].sum(axis=0)
        v = U["m"] / rho
        pl = (self.gam - 1.0) * U["ae"] / al - self.gam * self.pinf
        return {"al": al, "rl": rl, "v": v, "pl": pl, "A": U["A"], "chi": self.chi(al)}

    def derived(self, W, d=0):
        """Mixture quantities of a primitive state: density, stress, energies, wave speeds in direction d."""
        al, rl, v, pl, A, chi = W["al"], W["rl"], W["v"], W["pl"], W["A"], W["chi"]
        n = al.shape[1]
        rho = (al * rl).sum(axis=0)
        Y = al * rl / rho
        sig = -(al * pl).sum(axis=0)[:, None, None] * I3
        K = np.zeros((n, 3, 3))
        K[:, d, d] = (Y * self.gam * (pl + self.pinf) / rl).sum(axis=0)
        eel = np.zeros((self.N, n))
        for l, mat in enumerate(self.mats):
            if mat.mu > 0.0:
                g, g2, j2, ee, Sr = elastic(A[l], mat.mu, mat.rho0)
                sig += (chi * al[l] * rl[l])[:, None, None] * Sr
                eel[l] = chi * ee
                K += (chi * Y[l])[:, None, None] * acoustic_elastic(g, g2, j2, mat.mu, mat.rho0, d)
        t1, t2 = [i for i in range(3) if i != d]
        cl = np.sqrt(np.abs(K).sum(axis=2).max(axis=1))  # Gershgorin bound of the largest eigenvalue
        ct = np.sqrt(np.maximum(np.maximum(K[:, t1, t1], K[:, t2, t2]) + np.abs(K[:, t1, t2]), 0.0)) * chi
        E = (al * (pl + self.gam * self.pinf) / (self.gam - 1.0)).sum(axis=0) + (al * rl * eel).sum(axis=0) + 0.5 * rho * (v**2).sum(axis=0)
        return {"rho": rho, "sig": sig, "cl": cl, "ct": ct, "E": E, "eel": eel}

    def conservative(self, W):
        D = self.derived(W)
        return {
            "al": W["al"].copy(),
            "ar": W["al"] * W["rl"],
            "m": D["rho"] * W["v"],
            "E": D["E"],
            "ae": W["al"] * (W["pl"] + self.gam * self.pinf) / (self.gam - 1.0),
            "A": W["A"].copy(),
        }

    # ------------------------------------------------------------------------------------------------------------------
    # Riemann solvers at the faces (direction x). WL, WR: primitive states on both sides
    # ------------------------------------------------------------------------------------------------------------------
    def face(self, WL, WR):
        d = 0
        DL, DR = self.derived(WL, d), self.derived(WR, d)
        uL, uR = WL["v"][d], WR["v"][d]
        rL, rR = DL["rho"], DR["rho"]
        sL, sR = DL["sig"][:, d, d], DR["sig"][:, d, d]
        SL = np.minimum(uL - DL["cl"], uR - DR["cl"])
        SR = np.maximum(uL + DL["cl"], uR + DR["cl"])
        SM = (sL - sR + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / (rL * (SL - uL) - rR * (SR - uR))

        left = SM >= 0.0  # upwind side of the contact
        sup = np.where(left, SL >= 0.0, SR <= 0.0)  # supersonic: the upwind state itself

        def pick(a, b):
            return np.where(left, a, b)

        W = {k: pick(WL[k], WR[k]) if k != "A" else np.where(left[None, :, None, None], WL["A"], WR["A"]) for k in WL}
        rho, sig, E = pick(rL, rR), pick(sL, sR), pick(DL["E"], DR["E"])
        u = W["v"][d]
        S = pick(SL, SR)
        fac = np.where(sup, 1.0, (S - u) / (S - SM))
        ustar = np.where(sup, u, SM)

        # conservative part (longitudinal sub-system): F = F_K + S_K (U*_K - U_K)
        ar = W["al"] * W["rl"]
        F = {}
        F["ar"] = ar * fac * ustar
        mstar = rho * fac * W["v"]
        mstar[d] = rho * fac * ustar
        Fm = rho * u * W["v"]
        Fm[d] -= sig
        F["m"] = np.where(sup, Fm, Fm + S * (mstar - rho * W["v"]))
        Estar = fac * (E + (SM - u) * (rho * SM - sig / np.where(sup, 1.0, S - u)))
        F["E"] = np.where(sup, (E - sig) * u, (E - sig) * u + S * (Estar - E))

        # advected (non-conservative) quantities: flux u* phi*_K
        F["al"] = ustar * W["al"]
        rls = W["rl"] * fac
        pK = W["pl"] + self.pinf
        den = (self.gam - 1.0) * rls - (self.gam + 1.0) * W["rl"]
        hug = pK * ((self.gam - 1.0) * W["rl"] - (self.gam + 1.0) * rls) / np.where(den < 0.0, den, -1.0)
        pls = np.where((den < 0.0) & (hug > 0.0), hug, pK * fac**self.gam) - self.pinf  # Hugoniot, else isentrope
        F["ae"] = ustar * W["al"] * (pls + self.gam * self.pinf) / (self.gam - 1.0)
        As = W["A"].copy()
        As[..., d] *= fac[None, :, None]  # the normal column behaves like a density
        F["A"] = ustar[None, :, None, None] * As

        # shear sub-system: tangential tractions and velocities at the face
        ZL, ZR = rL * DL["ct"], rR * DR["ct"]
        Zs = ZL + ZR
        ok = Zs > 0.0
        Zs = np.where(ok, Zs, 1.0)
        vstar = np.zeros_like(WL["v"])
        vstar[d] = ustar
        for t in range(3):
            if t == d:
                continue
            tL, tR = DL["sig"][:, t, d], DR["sig"][:, t, d]
            vL, vR = WL["v"][t], WR["v"][t]
            vs = np.where(ok, (ZL * vL + ZR * vR + tR - tL) / Zs, 0.5 * (vL + vR))
            ts = np.where(ok, (ZR * tL + ZL * tR + ZL * ZR * (vR - vL)) / Zs, 0.0)
            F["m"][t] -= ts
            F["E"] -= ts * vs
            vstar[t] = vs
        return F, vstar

    # ------------------------------------------------------------------------------------------------------------------
    # space operator: returns dU/dt
    # ------------------------------------------------------------------------------------------------------------------
    def rhs(self, U, dx, order):
        ng = 2  # ghost cells on each side (transmissive: copies of the end cells)
        W = self.primitives(U)
        WL, WR = {}, {}
        for k, a in W.items():
            a = np.moveaxis(a, 1, -1) if k == "A" else a  # cell index last
            ag = np.concatenate([np.repeat(a[..., :1], ng, axis=-1), a, np.repeat(a[..., -1:], ng, axis=-1)], axis=-1)
            c = ag[..., 1:-1]  # cells -1 .. n (one ghost on each side)
            if order == 1 or k == "chi":
                s = np.zeros_like(c)
            else:
                dm = c - ag[..., :-2]
                dp = ag[..., 2:] - c
                s = np.where(dm * dp > 0.0, np.where(np.abs(dm) < np.abs(dp), dm, dp), 0.0)  # minmod
            lo, hi = (c + 0.5 * s)[..., :-1], (c - 0.5 * s)[..., 1:]  # left and right states of the n + 1 faces
            WL[k] = np.moveaxis(lo, -1, 1) if k == "A" else lo
            WR[k] = np.moveaxis(hi, -1, 1) if k == "A" else hi
        for Wf in (WL, WR):
            Wf["al"] = Wf["al"] / Wf["al"].sum(axis=0)
        F, vs = self.face(WL, WR)  # n + 1 faces

        d = 0
        dU = {}
        for k in ("ar", "m", "E", "al", "ae"):
            dU[k] = -(F[k][..., 1:] - F[k][..., :-1]) / dx
        dU["A"] = -(F["A"][:, 1:] - F["A"][:, :-1]) / dx
        du = (vs[:, 1:] - vs[:, :-1]) / dx  # (3, n): differences of the face velocities
        dU["al"] += W["al"] * du[d]
        dU["ae"] -= W["al"] * W["pl"] * du[d]
        A = W["A"]
        for t in range(3):
            if t == d:
                continue
            dU["A"][..., t] += A[..., t] * du[d][None, :, None]  # tangential columns: D/Dt = 0
            dU["A"][..., d] -= A[..., t] * du[t][None, :, None]  # normal column: + sum_t A[:, t] d(v_t)/dx
        return dU

    # ------------------------------------------------------------------------------------------------------------------
    # relaxations: fluidisation of the solid below alpha_lim, pressure relaxation, energy reset
    # ------------------------------------------------------------------------------------------------------------------
    def relax(self, U):
        al = U["al"] / U["al"].sum(axis=0)
        chi = self.chi(al)
        A = U["A"]
        for l in range(self.N):
            if self.solid[l]:
                J = np.abs(np.linalg.det(A[l])) ** (1.0 / 3.0)
                A[l] = np.where(chi[:, None, None] > 0.0, A[l], J[:, None, None] * I3)
        gam, pinf = self.gam, self.pinf
        p0 = (gam - 1.0) * U["ae"] / al - gam * pinf
        pmin = (-pinf).max() * (1.0 - 1e-12) + 1e-9
        p0 = np.maximum(p0, -pinf + 1e-9)  # phase pressures must stay above -pinf
        p = np.maximum((al * p0).sum(axis=0), pmin + 1e-6)
        for _ in range(60):  # Newton on f(p) = sum alpha_l tau_l(p)/tau_l0 - 1 (decreasing, convex)
            f = (al * (p0 + gam * pinf + (gam - 1.0) * p) / (gam * (p + pinf))).sum(axis=0) - 1.0
            df = -(al * (p0 + pinf) / (gam * (p + pinf) ** 2)).sum(axis=0)
            pn = p - f / df
            pn = np.where(pn <= pmin, 0.5 * (p + pmin), pn)
            if np.all(np.abs(pn - p) <= 1e-12 * (np.abs(p) + 1.0)):
                p = pn
                break
            p = pn
        al = al * (p0 + gam * pinf + (gam - 1.0) * p) / (gam * (p + pinf))
        al = al / al.sum(axis=0)
        U["al"] = al
        # reset of the phase energies from the total energy (mixture equation of state)
        rho = U["ar"].sum(axis=0)
        eel = np.zeros_like(al)
        chi = self.chi(al)
        for l, mat in enumerate(self.mats):
            if mat.mu > 0.0:
                eel[l] = chi * elastic(A[l], mat.mu, mat.rho0)[3]
        reh = U["E"] - 0.5 * (U["m"] ** 2).sum(axis=0) / rho - (U["ar"] * eel).sum(axis=0)
        pe = (reh - (al * gam * pinf / (gam - 1.0)).sum(axis=0)) / (al / (gam - 1.0)).sum(axis=0)
        if self.floor:
            # the energy-based pressure is kept only where admissible for every phase (p > -pinf_l); elsewhere (void
            # opening: gas close to vacuum, where tiny errors of the solid energy make it negative) the relaxed pressure
            # is kept and the total energy is made consistent with it (local loss of energy conservation)
            bad = pe <= pmin
            self.nfloor += int(bad.sum())
            p = np.where(bad, p, pe)
            U["ae"] = al * (p + gam * pinf) / (gam - 1.0)
            U["E"] = np.where(bad, U["ae"].sum(axis=0) + (U["ar"] * eel).sum(axis=0) + 0.5 * (U["m"] ** 2).sum(axis=0) / rho, U["E"])
        else:
            U["ae"] = al * (pe + gam * pinf) / (gam - 1.0)
        return U

    def step(self, U, dx, dt, order):
        def axpy(U, dU, a):
            return {k: U[k] + a * dU[k] for k in U}

        if order == 1:
            return self.relax(axpy(U, self.rhs(U, dx, 1), dt))
        U1 = self.relax(axpy(U, self.rhs(U, dx, 2), dt))
        U2 = axpy(U1, self.rhs(U1, dx, 2), dt)
        return self.relax({k: 0.5 * (U[k] + U2[k]) for k in U})

    def run(self, U, dx, tf, cfl=0.4, order=2, verbose=True):
        t, nt = 0.0, 0
        while t < tf * (1.0 - 1e-12):
            W = self.primitives(U)
            D = self.derived(W)
            dt = min(cfl * dx / np.max(np.abs(W["v"][0]) + D["cl"]), tf - t)
            U = self.step(U, dx, dt, order)
            t += dt
            nt += 1
            if not np.all(np.isfinite(U["E"])):
                raise FloatingPointError(f"NaN at step {nt}, t = {t}")
        if verbose:
            print(f"  {nt} steps, t = {t:.6e}")
        return U


# ----------------------------------------------------------------------------------------------------------------------
# cases
# ----------------------------------------------------------------------------------------------------------------------
def initial(model, x, alpha_solid, rho_l, v, p):
    """alpha_solid (n): volume fraction of phase 0; rho_l (N,n); v (3,n); p (n). Unstrained (A = I)."""
    n = x.size
    al = np.vstack([alpha_solid, 1.0 - alpha_solid])
    A = np.tile(I3, (model.N, n, 1, 1))
    W = {"al": al, "rl": rho_l, "v": v, "pl": np.vstack([p, p]), "A": A, "chi": model.chi(al)}
    return model.conservative(W)


def setup(case, n, eps_res=None):
    x = (np.arange(n) + 0.5) / n
    one = np.ones(n)
    v = np.zeros((3, n))
    p = 1e5 * one
    rho_s, rho_g = 8900.0 * one, 1.0 * one
    info = {}
    if case in ("advect", "slip-advect"):
        e = 1e-7 if eps_res is None else eps_res
        a = np.where(x < 0.5, 1.0 - e, e)
        v[0] = 1000.0
        if case == "slip-advect":
            v[1] = np.where(x < 0.5, 100.0, 0.0)
        tf = 0.15e-3
    elif case == "shear-slab":
        e = 1e-7 if eps_res is None else eps_res
        a = np.where((x > 0.15) & (x < 0.85), 1.0 - e, e)
        v[1] = np.where((x > 0.35) & (x < 0.65), 1000.0, 0.0)
        tf = 54.6e-6
        info["times"] = [54.6e-6, 180e-6]
    elif case == "solid-shear-tube":
        e = 1e-6 if eps_res is None else eps_res
        a = (1.0 - e) * one
        p = np.where(x < 0.5, 1e8, 1e5)
        v[1] = np.where(x < 0.5, 0.0, 100.0)
        tf = 0.46e-4
    elif case == "solid-gas-tube":
        e = 1e-6 if eps_res is None else eps_res
        a = np.where(x < 0.6, 1.0 - e, e)
        p = np.where(x < 0.6, 5e9, 1e5)
        rho_g = 50.0 * one
        tf = 87e-6
    elif case == "double-expansion":  # Favrie 2009 sec. 5.4: cavitation in copper
        e = 1e-4 if eps_res is None else eps_res
        a = (1.0 - e) * one
        v[0] = np.where(x < 0.5, -3000.0, 3000.0)
        tf = 88e-6
    elif case == "impact1d":  # copper slab at 800 m/s onto a copper slab at rest, in air
        e = 1e-4 if eps_res is None else eps_res
        a = np.where((x > 0.2) & (x < 0.6), 1.0 - e, e)
        v[0] = np.where((x > 0.2) & (x < 0.3), 800.0, 0.0)
        tf = 40e-6
    else:
        raise ValueError(case)
    return x, a, np.vstack([rho_s, rho_g]), v, p, tf, info


def dalembert_slab(x, t, ct, v0=1000.0, lo=0.15, hi=0.85, a=0.35, b=0.65):
    """Linear solution of the shear slab: even reflections at the free surfaces lo, hi. Returns v and the shear strain
    A[1][0] = -d(displacement_y)/dx."""
    L = hi - lo

    def box(xx):  # initial velocity extended evenly about lo and hi (period 2 L)
        s = np.mod(xx - lo, 2.0 * L)
        s = np.where(s > L, 2.0 * L - s, s) + lo
        return np.where((s > a) & (s < b), v0, 0.0)

    f, g = 0.5 * box(x - ct * t), 0.5 * box(x + ct * t)  # right- and left-going parts of the velocity
    inside = (x > lo) & (x < hi)
    v = np.where(inside, f + g, 0.0)
    a21 = np.where(inside, (f - g) / ct, 0.0)  # tau = Z (g - f), A[1][0] = -tau / mu
    return v, a21


def fields(model, U):
    W = model.primitives(U)
    D = model.derived(W)
    return {
        "rho": D["rho"], "alpha_s": W["al"][0], "u": W["v"][0], "v": W["v"][1], "p": (W["al"] * W["pl"]).sum(axis=0),
        "s11": D["sig"][:, 0, 0], "s12": D["sig"][:, 1, 0], "a21": U["A"][0][:, 1, 0] * W["chi"],
    }


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--cfl", type=float, default=0.4)
    ap.add_argument("--alim", type=float, default=0.5)
    ap.add_argument("--eps", type=float, default=None, help="residual volume fraction")
    ap.add_argument("--no-slip", action="store_true", help="disable the fluidisation (no slip treatment)")
    ap.add_argument("--no-floor", action="store_true", help="always reset the energies from the total energy")
    ap.add_argument("--tf", type=float, default=None)
    ap.add_argument("--plot", default=None)
    args = ap.parse_args(argv)

    model = Model([COPPER, AIR], alim=args.alim, slip=not args.no_slip, floor=not args.no_floor)
    x, a, rl, v, p, tf, info = setup(args.case, args.n, args.eps)
    if args.tf is not None:
        tf = args.tf
    U = initial(model, x, a, rl, v, p)
    U = model.run(U, 1.0 / args.n, tf, args.cfl, args.order)
    f = fields(model, U)
    for k, val in f.items():
        print(f"  {k:8s} min {val.min(): .6e} max {val.max(): .6e}")
    print(f"  cells x stages where the relaxed pressure was kept: {model.nfloor}")
    if args.case == "shear-slab":
        ct = np.sqrt(COPPER.mu / COPPER.rho0)
        vex, aex = dalembert_slab(x, tf, ct)
        print(f"  L1 error on v vs linear solution: {np.abs(f['v'] - vex).mean():.4e} m/s;"
              f" max |v| in the air: {np.abs(f['v'][(x < 0.14) | (x > 0.86)]).max():.4e}")
    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        keys = list(f)
        fig, axs = plt.subplots(2, 4, figsize=(16, 7))
        for ax, k in zip(axs.ravel(), keys):
            ax.plot(x, f[k], "k.-", ms=2, lw=0.6)
            if args.case == "shear-slab" and k in ("v", "a21"):
                ax.plot(x, vex if k == "v" else aex, "r-", lw=0.8)
            ax.set_title(k)
            ax.grid(alpha=0.3)
        fig.suptitle(f"{args.case}  n={args.n} order={args.order} slip={'on' if model.slip else 'off'} t={tf:.3e}")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=110)
        print("  wrote", args.plot)
    return x, f


if __name__ == "__main__":
    main()
