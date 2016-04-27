# running the Analytical solution for SAG
using Roots, Plots, ProgressMeter, DataFrames, Dierckx, JLD
include("frac_flow_funcs.jl")
(xt_shock, sw_shock, xt_prf, sw_prf, t, p_inj)=frac_flow_sag(
fmmob=2.5e4, epdry=1.0e5, fmdry=0.29,
muw=1e-3, mug=0.02e-3, ut=1e-5, phi=0.2,
k=1e-12, swc=0.2, sgr=0.2, krg0=0.94, ng=1.3, krw0=0.2,
  nw=4.2, sw0=1.0, sw_inj=0.0, L=1, pv_inj=5)
plot(xt_prf, sw_prf)
(krw, krg, dkrwdsw, dkrgdsw, krgf, dkrgfdsw)=corey_rel_perms(fmmob=2.5e4, epdry=1.0e5, fmdry=0.29, swc=0.2, sgr=0.2, krg0=0.94, ng=1.3, krw0=0.2,
  nw=4.2)
plot(sw_prf, krw(sw_prf)/1e-3+krgf(sw_prf)/0.02e-3)
