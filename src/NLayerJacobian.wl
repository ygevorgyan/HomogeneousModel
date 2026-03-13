(* ::Package:: *)
(* ================================================================= *)
(* NLayerJacobian.wl — Unified 6x6 propagation for N-layer bodies   *)
(* ================================================================= *)
(* Module 5 of the modular k2(s) pipeline.                           *)
(*                                                                     *)
(* Formulation (unified 6x6):                                         *)
(*   All layers (2..N) propagated through full 6x6 Y-matrices.       *)
(*   Liquid layers use muHat ~ 10^-20 regularization (no special     *)
(*   Clairaut propagator needed).  Works for arbitrary layer configs: *)
(*   any combination of solid/liquid layers in any position.          *)
(*   W1    = 3x3 secular matrix (3 surface BCs x 3 IC amplitudes)    *)
(*   det(W1) = 0 gives secular poles                                  *)
(*   k2(s) via LinearSolve on W1                                     *)
(*                                                                     *)
(* Exports:                                                           *)
(*   nLayerShellProp[rk,rhok,muk,etak,gk,s]       -> 6x6             *)
(*   nLayerW1[rk,rhok,muk,etak,gk,s]              -> {W1,Gam}        *)
(*   nLayerDetNum[s, rk,rhok,muk,etak]             -> scalar         *)
(*   nLayerK2Num[s, rk,rhok,muk,etak]              -> scalar         *)
(*   nLayerResidues[poles, rk,rhok,muk,etak]       -> {r1,...}       *)
(*                                                                     *)
(* Depends on: LayeredBodyConstants`, YMatrix`                        *)
(* ================================================================= *)

BeginPackage["NLayerJacobian`",
  {"LayeredBodyConstants`", "YMatrix`"}];

$NLwp::usage =
  "$NLwp is the working precision for all numerical evaluations \
(default 50). Needed to handle ill-conditioned propagator matrices \
and FindRoot convergence.";

nLayerShellProp::usage =
  "nLayerShellProp[rk, rhok, muk, etak, gk, sval] returns the 6x6 \
shell propagator Bshell = P_N . ... . P_2 for all layers above the \
inner core (j = 2 to N). Liquid layers use muHat regularization.";

nLayerW1::usage =
  "nLayerW1[rk, rhok, muk, etak, gk, sval] returns {W1, Gam} \
where W1 is the 3x3 secular matrix and Gam is 6x3. \
The secular equation is Det[W1] = 0.";

nLayerDetNum::usage =
  "nLayerDetNum[sval, rk, rhok, muk, etak] returns the numerical \
determinant Det[W1(s)]. Gravity is computed internally.";

nLayerK2Num::usage =
  "nLayerK2Num[sval, rk, rhok, muk, etak] returns k2(s) via \
LinearSolve on the 3x3 secular matrix W1.";

nLayerResidues::usage =
  "nLayerResidues[poles, rk, rhok, muk, etak] returns the numerical \
residues of k2(s) at the given poles via the limit formula \
r_j = delta * k2(s_j + delta).";

Begin["`Private`"];

$NLwp = 50;

(* ---- Internal helpers ------------------------------------------------ *)

(* Maxwell muhat at a specific numeric s-value, with liquid regularization *)
muHatNum[muk_, etak_, sval_] := Table[
  If[muk[[j]] < 1,
    SetPrecision[10^-20, $NLwp],
    muk[[j]] sval / (sval + muk[[j]] / etak[[j]])],
  {j, Length[muk]}];

(* Promote all parameter arrays to working precision *)
promote[rk_, rhok_, muk_, etak_, gk_, sval_] :=
  {SetPrecision[rk, $NLwp], SetPrecision[rhok, $NLwp],
   SetPrecision[muk, $NLwp], SetPrecision[etak, $NLwp],
   SetPrecision[gk, $NLwp], SetPrecision[sval, $NLwp]};

(* ---- Constant matrices (module-level) -------------------------------- *)

(* P1: extraction {R,S,Q} from 6-component state vector (3x6) *)
$P1 = {{0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 1}};

(* P5: extraction {Phi} from 6-component state vector (1x6) *)
$P5 = {{0, 0, 0, 0, 1, 0}};

(* ---- Shell propagator (all layers j=2..N) ----------------------------- *)

nLayerShellProp[rk0_, rhok0_, muk0_, etak0_, gk0_, sval0_] :=
  Module[{rk, rhok, muk, etak, gk, sval, nn, muhats, Bsh, YIprev},
    {rk, rhok, muk, etak, gk, sval} =
        promote[rk0, rhok0, muk0, etak0, gk0, sval0];
    nn = Length[rk];
    muhats = muHatNum[muk, etak, sval];
    Bsh = IdentityMatrix[6];
    Do[
      YIprev = LinearSolve[
        Ymat[rk[[j - 1]], rhok[[j]], muhats[[j]], gk[[j - 1]]],
        IdentityMatrix[6]];
      Bsh = Ymat[rk[[j]], rhok[[j]], muhats[[j]], gk[[j]]] .
        YIprev . Bsh,
      {j, 2, nn}];
    Bsh
  ];

(* ---- W1 (3x3 secular matrix) ---------------------------------------- *)

nLayerW1[rk0_, rhok0_, muk0_, etak0_, gk0_, sval0_] :=
  Module[{rk, rhok, muk, etak, gk, sval, muhats,
          IcR1, Bsh, Gam, W1},
    {rk, rhok, muk, etak, gk, sval} =
        promote[rk0, rhok0, muk0, etak0, gk0, sval0];
    muhats = muHatNum[muk, etak, sval];

    (* Inner core basis at r1 (6x3) *)
    IcR1 = IcMat[rk[[1]], rhok[[1]], muhats[[1]], gk[[1]]];

    (* Shell propagator for all layers j=2..N *)
    Bsh = nLayerShellProp[rk, rhok, muk, etak, gk, sval];

    (* Propagated basis at surface (6x3) *)
    Gam = Bsh . IcR1;

    (* Surface BCs: W1 = P1 . Gam (3x3) *)
    W1 = $P1 . Gam;

    {W1, Gam}
  ];

(* ---- Det[W1] and k2(s) ---------------------------------------------- *)

nLayerDetNum[sval_?NumericQ, rk_, rhok_, muk_, etak_] :=
  Module[{gk, result},
    gk = computeGravity[rk, rhok];
    result = nLayerW1[rk, rhok, muk, etak, gk, sval];
    Det[result[[1]]]
  ];

nLayerK2Num[sval_?NumericQ, rk_, rhok_, muk_, etak_] :=
  Module[{gk, nn, result, W1, Gam, rhsVec, sol, PhiA},
    gk = computeGravity[rk, rhok];
    nn = Length[rk];
    result = nLayerW1[rk, rhok, muk, etak, gk, sval];
    W1 = result[[1]];
    Gam = result[[2]];

    (* RHS: {0, 0, (2l+1)/R} where R = rk[[nn]] *)
    rhsVec = SetPrecision[{0, 0, (2 $lDegree + 1)/rk[[nn]]}, $NLwp];

    (* Solve W1 . sol = rhsVec *)
    sol = LinearSolve[W1, rhsVec];

    (* k2 = (P5 . Gam . sol)[[1]] - 1 *)
    PhiA = $P5 . Gam . sol;
    PhiA[[1]] - 1
  ];

(* Residues via limit: r_j = delta * k2(s_j + delta) *)
nLayerResidues[poles_, rk_, rhok_, muk_, etak_] :=
  Table[
    Module[{delta},
      delta = SetPrecision[10^-6 Abs[poles[[j]]], $NLwp];
      Re[delta nLayerK2Num[poles[[j]] + delta, rk, rhok, muk, etak]]],
    {j, Length[poles]}];

End[];

EndPackage[];
