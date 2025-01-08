using DrWatson
@quickactivate "KLAP"

using LinearAlgebra, ControlSystemsBase, RobustAndOptimalControl
using DataFrames, Plots, LaTeXStrings, JuMP, Optim
using Passivation

import Passivation: ispassive

# Smartphone model
Σ = wload(datadir("smartphone.jld2"))["Σ"];
Σref = wload(datadir("smartphone_ref.jld2"))["Σ"];
h2_error_ref = h2norm(Σ - Σref)

# Diagonalize system
F = eigen(Σ.A);
Σdiag = ss(diagm(F.values), F.vectors \ Σ.B, Σ.C * F.vectors, Σ.D);
# h2norm(Σ - Σdiag)

# ΔD = 8.0 is necessary to make the computation of the ARE feasible
L0, ΔD = klap_inital_guess(Σdiag, 8.0);
Σp, res = passivate(Σdiag, :klap, L0; verbose=true, show_every=100);

h2_error = h2norm(Σ - Σp)

@info "Smartphone example: solve time: $res.time_run,\t iterations: $res.iterations,\t H2-error: $h2_error"

improvement = (h2_error_ref - h2_error) / h2_error_ref
@info "Relative improvement: $(improvement*100)%"

p1 = sigmaplot([Σ, Σref, Σp]; extrema=true, label=["Non-passive" "Reference" "KLAP"], lc=[:black 1 2], ls=[:solid :dash :dash], xlims=[1e6,1e12]);
p2 = sigmaplot([Σ - Σref, Σ - Σp]; extrema=true, label=["Reference" "KLAP"], lc=[1 2], ls=[:dash :dash], xlims=[1e6,1e12], title="Sigma Error Plot");

plot(p1, p2, layout=(2,1))