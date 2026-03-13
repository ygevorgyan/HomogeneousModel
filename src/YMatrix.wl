(* ::Package:: *)
(* ================================================================= *)
(* YMatrix.wl — 6x6 fundamental matrix Y and 6x3 inner core Ic      *)
(* ================================================================= *)
(* Module 2 of the modular k2(s) pipeline.                           *)
(* Sabadini Eq. 2.42 specialized to l=2, incompressible.             *)
(* State vector ordering: {U, V, R, S, Phi, Q}.                     *)
(*                                                                     *)
(* Exports: Ymat[r, rho, muhat, g]  -> 6x6                          *)
(*          IcMat[r, rho, muhat, g] -> 6x3                          *)
(* Depends on: LayeredBodyConstants` (for $GravConst)                *)
(* ================================================================= *)

BeginPackage["YMatrix`", {"LayeredBodyConstants`"}];

Ymat::usage =
  "Ymat[r, rho, muhat, g] returns the 6x6 fundamental matrix Y \
(Sabadini Eq. 2.42, l=2, incompressible). Arguments: r = radius, \
rho = density, muhat = complex shear modulus, g = gravity.";

IcMat::usage =
  "IcMat[r, rho, muhat, g] returns the 6x3 inner core matrix \
(first 3 columns of Y: solutions regular at r=0). \
Eq. A.3 of Matsuyama et al. (2018).";

Begin["`Private`"];

Ymat[r_, rho_, muhat_, g_] := {
  {r^3/7, r, 0, 1/(2 r^2), 1/r^4, 0},
  {5 r^3/42, r/2, 0, 0, -1/(3 r^4), 0},
  {r^2 (g rho r - muhat)/7, g rho r + 2 muhat, rho r^2,
   (g rho r - 6 muhat)/(2 r^3), (g rho r - 8 muhat)/r^5, rho/r^3},
  {8 muhat r^2/21, muhat, 0, muhat/(2 r^3), 8 muhat/(3 r^5), 0},
  {0, 0, r^2, 0, 0, 1/r^3},
  {4 Pi $GravConst rho r^3/7, 4 Pi $GravConst rho r, 5 r,
   2 Pi $GravConst rho/r^2, 4 Pi $GravConst rho/r^4, 0}
};

IcMat[r_, rho_, muhat_, g_] := {
  {r^3/7, r, 0},
  {5 r^3/42, r/2, 0},
  {r^2 (g rho r - muhat)/7, g rho r + 2 muhat, rho r^2},
  {8 muhat r^2/21, muhat, 0},
  {0, 0, r^2},
  {4 Pi $GravConst rho r^3/7, 4 Pi $GravConst rho r, 5 r}
};

End[];

EndPackage[];
