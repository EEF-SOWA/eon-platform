// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>
#include <algorithm>

namespace eon {

// ---------------------------------------------------------------------------
// Koren triode model + REAL Newton-Raphson plate-voltage solve.
//
// Koren model (12AX7-flavoured defaults):
//   E1 = (Vp/Kp) * ln(1 + exp(Kp*(1/mu + Vg/sqrt(Kvb + Vp^2))))
//   Ip = E1^ex / Kg                                        (E1 > 0)
//
// Circuit: plate supply B+ through load Rload.
//   f(Vp) = (B+ - Vp)/Rload - Ip(Vg, Vp) = 0
//
// Newton step uses the ANALYTIC dIp/dVp (the original used a fabricated
// derivative `Ip*0.01` and solved an equation that wasn't the circuit).
// Solution is warm-started from the previous sample — converges in ~2 iters.
// ---------------------------------------------------------------------------

struct TriodeStage
{
    // Koren parameters
    float mu  = 100.f;   // amplification factor
    float ex  = 1.4f;    // exponent
    float kg  = 1060.f;  // (roughly) inverse transconductance scale
    float kp  = 600.f;
    float kvb = 300.f;

    // Circuit
    float Bplus  = 250.f;     // HT supply [V]
    float Rload  = 100000.f;  // total plate+cathode load [ohm]
    float biasVg = -1.5f;     // DC grid bias [V]
    float inVolts = 3.0f;     // input sensitivity: full-scale -> volts of Vg swing

    int iterations = 2;

    // --- optional self-biasing cathode (Rk || Ck to ground) ---------------
    // Off by default: cathodeRk == 0 keeps cathodeVk == 0 and the output is
    // bit-identical to the grounded-cathode version.
    //
    // Enabled via setCathode(): the cathode rides on Rk||Ck so the average
    // plate current lifts Vk, which reduces Vgk and self-biases the tube
    // colder — program-dependent "cathode sag" compression (ported from the
    // Tube comp TriodeStage, which uses the same mechanism for vari-mu feel).
    float  cathodeRk   = 0.f;     // cathode resistor [ohm], 0 = grounded
    float  cathodeCk   = 25e-6f;  // bypass capacitor [F]
    double cathodePole = 0.0;     // exp(-1/(Rk*Ck*rate)) — set by setCathode
    double cathodeVk   = 0.0;     // tracked cathode voltage [V]
    float  lastIp      = 0.f;     // plate current from last solve [A]

    // state
    float VpPrev = 150.f;
    float Vq = 150.f;         // quiescent plate voltage (solved at reset)
    float outScale = 1.f;

    // --- optional output coupling cap (one-pole high-pass, off by default) ---
    // A real stage couples into the next one through a cap, so large swings
    // that shift the plate's DC level decay instead of persisting. fcHz <= 0
    // disables it (output unchanged). rateHz = processing rate.
    double coupleR = 0.0, coupleX1 = 0.0, coupleY1 = 0.0;
    void setOutputCoupling (double fcHz, double rateHz)
    {
        coupleR = (fcHz > 0.0 && rateHz > 0.0)
            ? std::exp (-2.0 * 3.14159265358979323846 * fcHz / rateHz) : 0.0;
        coupleX1 = coupleY1 = 0.0;
    }

    // Fit the Rk||Ck cathode network. rateHz is the *processing* rate —
    // under oversampling pass the oversampled rate, not the host rate.
    // Call before reset() so the quiescent solve sees the fitted network.
    void setCathode (float rkOhms, float ckFarads, double rateHz)
    {
        cathodeRk = rkOhms;
        cathodeCk = ckFarads;
        cathodePole = (rkOhms > 0.f && ckFarads > 0.f && rateHz > 0.0)
            ? std::exp (-1.0 / ((double) rkOhms * ckFarads * rateHz)) : 0.0;
    }

    void reset()
    {
        const int processingIterations = iterations;
        iterations = 16;
        VpPrev = 150.f;
        cathodeVk = 0.0;
        lastIp = 0.f;
        coupleX1 = coupleY1 = 0.0;

        // Reset can spend more work than the per-sample solver to find a
        // consistent DC operating point without changing the processing budget.
        const bool selfBiased = cathodeRk > 0.f && cathodePole > 0.0;
        double Vk = 0.0;
        if (selfBiased)
        {
            constexpr int maxBisectionIterations = 48;
            constexpr double cathodeTolerance = 1e-5;
            const auto residualAt = [this](double candidateVk)
            {
                // Give every residual evaluation the same plate-solver seed;
                // the final selected point is solved again below.
                VpPrev = 150.f;
                solveVp (biasVg, (float) candidateVk);
                return (double) lastIp * cathodeRk - candidateVk;
            };

            double lowerVk = 0.0;
            double upperVk = std::max (0.0, (double) Bplus);
            double lowerResidual = residualAt (lowerVk);
            double upperResidual = residualAt (upperVk);
            double selectedVk = std::abs (lowerResidual) <= std::abs (upperResidual)
                ? lowerVk : upperVk;
            double selectedResidual = std::min (std::abs (lowerResidual),
                                                std::abs (upperResidual));

            const bool bracketed = (lowerResidual <= 0.0 && upperResidual >= 0.0)
                                || (lowerResidual >= 0.0 && upperResidual <= 0.0);
            if (bracketed && selectedResidual > cathodeTolerance)
            {
                for (int i = 0; i < maxBisectionIterations; ++i)
                {
                    const double midpoint = 0.5 * (lowerVk + upperVk);
                    if (midpoint == lowerVk || midpoint == upperVk)
                        break;

                    const double midpointResidual = residualAt (midpoint);
                    if (std::abs (midpointResidual) < selectedResidual)
                    {
                        selectedVk = midpoint;
                        selectedResidual = std::abs (midpointResidual);
                    }
                    if (selectedResidual <= cathodeTolerance)
                        break;

                    if ((lowerResidual <= 0.0 && midpointResidual <= 0.0)
                        || (lowerResidual >= 0.0 && midpointResidual >= 0.0))
                    {
                        lowerVk = midpoint;
                        lowerResidual = midpointResidual;
                    }
                    else
                    {
                        upperVk = midpoint;
                        upperResidual = midpointResidual;
                    }
                }
            }
            // If the physical interval does not bracket a root, keep the
            // endpoint with the smaller residual instead of treating an
            // unbracketed iteration as converged. At the iteration cap,
            // selectedVk retains the lowest residual sampled so far.
            Vk = selectedVk;
        }

        cathodeVk = Vk;
        VpPrev = 150.f;
        Vq = solveVp (biasVg, (float) cathodeVk);
        outScale = 1.f / std::max (40.f, Bplus - Vq);
        VpPrev = Vq;
        iterations = processingIterations;
    }

    // Templated so the per-sample solver can run in double while callers that
    // pass floats keep the float arithmetic they always had.
    template <typename T>
    inline T korenIp (T Vg, T Vp, T* dIp_dVp = nullptr) const
    {
        const T Kp = (T) kp, Ex = (T) ex, Kg = (T) kg;
        const T S  = std::sqrt ((T) kvb + Vp * Vp);
        const T u  = Kp * (T (1) / (T) mu + Vg / S);
        // softplus(u) = ln(1+e^u), numerically stable
        const T sp = u > T (20) ? u : std::log1p (std::exp (u));
        const T E1 = Vp / Kp * sp;
        const T Ip = E1 > T (0) ? std::pow (E1, Ex) / Kg : T (0);

        if (dIp_dVp != nullptr)
        {
            const T sig = u > T (20) ? T (1) : T (1) / (T (1) + std::exp (-u));
            // dE1/dVp = sp/kp + (Vp/kp)*sig*(du/dVp), du/dVp = -kp*Vg*Vp/S^3
            const T dE1 = sp / Kp - Vg * Vp * Vp * sig / (S * S * S);
            *dIp_dVp = E1 > T (0) ? Ex * std::pow (E1, Ex - T (1)) / Kg * dE1 : T (0);
        }
        return Ip;
    }

    // Solve f(Vp) = (B+ - Vp)/Rload - Ip(Vg - Vk, Vp - Vk) = 0.
    // Vk is the (slowly varying) cathode voltage; the dVpk/dVp term is 1 so
    // the analytic derivative applies unchanged. Vk == 0 -> legacy solve.
    inline float solveVp (float Vg, float Vk = 0.f)
    {
        // The Newton loop runs in double; state and results stay float.
        const double vg = (double) Vg - (double) Vk, vk = Vk;
        const double hi = Bplus, rl = Rload;
        double Vp = VpPrev;
        for (int i = 0; i < iterations; ++i)
        {
            double dIp;
            const double Ip = korenIp<double> (vg, Vp - vk, &dIp);
            const double f  = (hi - Vp) / rl - Ip;
            const double df = -1.0 / rl - dIp;     // strictly negative -> safe
            Vp = std::clamp (Vp - f / df, 0.0, hi);
        }
        lastIp = (float) korenIp<double> (vg, Vp - vk);
        return (float) Vp;
    }

    // Returns normalized plate swing: >0 for positive grid drive.
    inline float process (float x)
    {
        const float Vg = biasVg + x * inVolts;
        const float Vp = solveVp (Vg, (float) cathodeVk);
        VpPrev = Vp;
        // Rk||Ck one-pole tracks Ip*Rk -> self-bias sag (off when unconfigured)
        if (cathodeRk > 0.f && cathodePole > 0.0)
            cathodeVk += (1.0 - cathodePole) * ((double) lastIp * cathodeRk - cathodeVk);
        const float y = (Vq - Vp) * outScale;
        if (coupleR <= 0.0)
            return y;
        const double yc = (double) y - coupleX1 + coupleR * coupleY1;
        coupleX1 = y; coupleY1 = yc;
        return (float) yc;
    }
};

} // namespace eon
