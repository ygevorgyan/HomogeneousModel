(* ::Package:: *)
(* ================================================================= *)
(* LayeredBodyConstants.wl — Physical constants and generic helpers   *)
(*   for an N-layer viscoelastic body                                *)
(* ================================================================= *)
(* Module 1 of the modular k2(s) pipeline.                           *)
(* Exports: $GravConst, $lDegree, fmtSci[],                          *)
(*   computeGravity[], muHat[], computeMoI[],                        *)
(*   computePrefactor[], computeFconv[]                              *)
(* ================================================================= *)

BeginPackage["LayeredBodyConstants`"];

$GravConst::usage =
  "$GravConst is the gravitational constant G = 6674/10^14 m^3 kg^-1 s^-2.";
$lDegree::usage =
  "$lDegree is the spherical harmonic degree (2).";

computeGravity::usage =
  "computeGravity[rk, rhok] computes gravity g(r) = G*M(r)/r^2 at each \
layer interface. Returns a list of length Length[rk].";

muHat::usage =
  "muHat[mu, eta, s] returns the scalar complex rigidity for a Maxwell \
material: mu*eta*s/(mu + eta*s). All arguments are symbolic.";

computeMoI::usage =
  "computeMoI[rk, rhok] computes the mean moment of inertia \
I = (8*Pi/15) * Sum[rho_j*(r_j^5 - r_{j-1}^5)].";

computePrefactor::usage =
  "computePrefactor[rk, rhok] computes the mapping prefactor \
3*I*G/R^5 (Gevorgyan et al. 2023, Eq. 2).";

computeFconv::usage =
  "computeFconv[rk, rhok] computes the unit conversion factor \
F = 5*rhobar*R^2/38 from compliance space to Pa (for l=2).";

fmtSci::usage =
  "fmtSci[x, ndig] formats a number x in scientific notation with \
ndig significant digits, returning a string like \"1.23e4\".";

Begin["`Private`"];

$GravConst = 6674/10^14;
$lDegree = 2;

computeGravity[rk_, rhok_] := Module[{gka, nn},
  nn = Length[rk];
  gka = 4 Pi Flatten[{
    rhok[[1]] rk[[1]]^3/3,
    Table[
      rhok[[1]] rk[[1]]^3/3 +
      Sum[rhok[[j]] (rk[[j]]^3/3 - rk[[j-1]]^3/3), {j, 2, k}],
      {k, 2, nn}]
  }];
  $GravConst Table[gka[[j]]/rk[[j]]^2, {j, 1, nn}]
];

muHat[mu_, eta_, s_] := mu eta s/(mu + eta s);

computeMoI[rk_, rhok_] := Module[{rkExt, nn},
  nn = Length[rk];
  rkExt = Prepend[rk, 0];
  (8 Pi/15) Sum[rhok[[j]] (rk[[j]]^5 - rkExt[[j]]^5), {j, 1, nn}]
];

computePrefactor[rk_, rhok_] := Module[{Imoi, R},
  Imoi = computeMoI[rk, rhok];
  R = rk[[-1]];
  3 Imoi $GravConst/R^5
];

computeFconv[rk_, rhok_] := Module[{rkExt, nn, Mtotal, R, rhobar},
  nn = Length[rk];
  rkExt = Prepend[rk, 0];
  Mtotal = (4 Pi/3) Sum[rhok[[j]] (rk[[j]]^3 - rkExt[[j]]^3), {j, 1, nn}];
  R = rk[[-1]];
  rhobar = 3 Mtotal/(4 Pi R^3);
  5 rhobar R^2/38
];

fmtSci[x_, ndig_] := Module[{xr, e, m},
  xr = Re[N[x]];
  If[!NumberQ[xr] || xr == 0, Return["0"]];
  If[Abs[xr] === Infinity || !NumericQ[Abs[xr]], Return[ToString[xr]]];
  e = Floor[Log10[Abs[xr]]];
  m = xr/10.^e;
  ToString[NumberForm[m, {ndig, ndig - 1}]] <> "e" <> ToString[e]
];

End[];

EndPackage[];
