(* ::Package:: *)
(* ================================================================= *)
(* HomogeneousRheology.wl — Equivalent homogeneous body from layers  *)
(* ================================================================= *)
(* Computes the equivalent homogeneous Voigt rheology from a layered *)
(* model: Love numbers, secular poles, residues, and generalized     *)
(* Voigt element parameters (gamma, alpha, eta, {tau_i, alpha_i,     *)
(* eta_i}).                                                           *)
(*                                                                     *)
(* Extracted from run_pipeline.wl Section 1 (computeVoigt).          *)
(*                                                                     *)
(* Exports: computeHomogeneous[rk, rhok, muk, etak, omega, opts]     *)
(*                                                                     *)
(* Depends on: LayeredBodyConstants` (computePrefactor, computeFconv),*)
(*             NLayerJacobian` (nLayerK2Num, nLayerDetNum,           *)
(*                              nLayerResidues)                      *)
(* ================================================================= *)

BeginPackage["HomogeneousRheology`", {"LayeredBodyConstants`", "NLayerJacobian`"}];

computeHomogeneous::usage =
  "computeHomogeneous[rk, rhok, muk, etak, omega, opts] computes the \
equivalent homogeneous Voigt rheology from a layered model. Returns an \
Association with keys: kE, kf, k2Re, k2Im, gamma, alpha, eta, \
prefactor, Fconv, nPoles, poles, residues, nVoigt, tauV, alphaV, etaV. \
Options: \"WorkingPrecision\" (default 50), \"GridLogMin\" (default -17), \
\"GridLogMax\" (default -5), \"GridPoints\" (default 500).";

Begin["`Private`"];

Options[computeHomogeneous] = {
  "WorkingPrecision" -> 50,
  "GridLogMin" -> -17,
  "GridLogMax" -> -5,
  "GridPoints" -> 500
};

computeHomogeneous[rk_, rhok_, muk_, etak_, omega_, opts:OptionsPattern[]] :=
  Module[
  {wp, gridLogMin, gridLogMax, gridN,
   kEval, kfval, k2val, k2Re, k2Im, Aj,
   prefactor, Fconv,
   gammac, alphac, gammaPa, alphaPa, eta0c, eta0Pa,
   sGrid, detGrid, brackets, poles9, nPoles, residj,
   Jres0val, nVoigt, etaVc, alphaVc, alphaVPa, etaVPa, tauV,
   numPoly, denPoly, JnumPoly, JdenFullPoly, JdenRPoly, ss,
   Jpoles, Jresidues, JdenRPrime},

  (* Unpack options *)
  wp = OptionValue["WorkingPrecision"];
  gridLogMin = OptionValue["GridLogMin"];
  gridLogMax = OptionValue["GridLogMax"];
  gridN = OptionValue["GridPoints"];

  (* 1. kE: elastic Love number at large s *)
  kEval = Re[nLayerK2Num[SetPrecision[10^10, wp], rk, rhok, muk, etak]];

  (* 2. k2(iw) at tidal frequency *)
  k2val = nLayerK2Num[SetPrecision[I omega, wp], rk, rhok, muk, etak];
  k2Re = Re[k2val];
  k2Im = Im[k2val];

  (* 3. Find all secular poles via grid search + Brent refinement.
        Expected number of modes for l=2, incompressible
        (Sabadini, Vermeersen & Cambiotti 2016, Sect. 1.8):
          - Each solid-solid interface: up to 3 modes
              1 buoyancy (if drho =/= 0) + 2 transient (if Maxwell times differ)
          - Each liquid-solid interface: 1 buoyancy mode
          - Free surface above viscoelastic layer: 1 mode
        So for generic parameters: M = 3*Nss + Nls + 1
        where Nss = number of solid-solid internal interfaces,
              Nls = number of liquid-solid internal interfaces.
        Verified: Moon=9 (2 s-s + 2 l-s + surface = 6+2+1),
                  Europa=4 (0 s-s + 3 l-s + surface = 0+3+1),
                  Io=6 (2 s-s + 1 l-s + surface, but drho=0 at s-s -> 4+1+1),
                  Enceladus=3 (0 s-s + 2 l-s + surface = 0+2+1). *)
  sGrid = -10^Subdivide[SetPrecision[gridLogMin, wp], SetPrecision[gridLogMax, wp], gridN];
  detGrid = Table[
    Re[nLayerDetNum[sv, rk, rhok, muk, etak]],
    {sv, sGrid}];

  (* Sign changes bracket roots *)
  brackets = {};
  Do[
    If[detGrid[[i]] detGrid[[i + 1]] < 0,
      AppendTo[brackets, {sGrid[[i]], sGrid[[i + 1]]}]],
    {i, Length[detGrid] - 1}];

  (* Refine each bracket with FindRoot/Brent *)
  poles9 = Table[
    ss /. Quiet[FindRoot[
      Re[nLayerDetNum[ss, rk, rhok, muk, etak]],
      {ss, br[[1]], br[[2]]},
      WorkingPrecision -> wp, Method -> "Brent"], {FindRoot::brdig}],
    {br, brackets}];
  poles9 = Sort[poles9, Abs[#1] > Abs[#2] &];
  nPoles = Length[poles9];

  (* 4. Residues via Module 5 limit formula *)
  residj = nLayerResidues[poles9, rk, rhok, muk, etak];

  (* 5. kf from partial fractions (avoids s~0 singularity) *)
  Aj = Table[-residj[[j]]/poles9[[j]], {j, nPoles}];
  kfval = kEval + Total[Aj];

  (* 6. Mapping constants from Module 1 *)
  prefactor = computePrefactor[rk, rhok];
  Fconv = computeFconv[rk, rhok];

  gammac = prefactor/kfval;
  alphac = prefactor/kEval - gammac;
  gammaPa = gammac Fconv;
  alphaPa = alphac Fconv;

  (* 7. Voigt inversion via rational function reconstruction *)
  (* k2(s) = kE + Sum[rj/(s-sj)] rewritten as numPoly/denPoly *)
  numPoly = kEval Product[ss - poles9[[j]], {j, nPoles}] +
    Sum[residj[[j]] Product[If[k == j, 1, ss - poles9[[k]]], {k, nPoles}],
      {j, nPoles}];
  denPoly = Product[ss - poles9[[j]], {j, nPoles}];

  (* J(s) = k2/(prefactor - gammac*k2) => numerator=numPoly,
     denominator = prefactor*denPoly - gammac*numPoly *)
  JdenFullPoly = Expand[prefactor denPoly - gammac numPoly];

  (* Factor out the ss=0 root (guaranteed since gammac = prefactor/kf) *)
  JdenRPoly = Sum[
    Coefficient[JdenFullPoly, ss, k] ss^(k - 1),
    {k, 1, Exponent[JdenFullPoly, ss]}];

  JnumPoly = Expand[numPoly];

  (* Residue of J at ss=0 => isolated dashpot: eta0 = 1/Jres0 *)
  Jres0val = (JnumPoly /. ss -> 0)/(JdenRPoly /. ss -> 0);
  eta0c = 1/Jres0val;
  eta0Pa = eta0c Fconv;

  (* Find Voigt element poles: roots of JdenRPoly *)
  (* Boost precision to 2*wp to handle polynomial cancellation *)
  Jpoles = ss /. Quiet[NSolve[SetPrecision[JdenRPoly, 2 wp] == 0, ss,
    WorkingPrecision -> 2 wp], {NSolve::precw}];
  Jpoles = Sort[Re[Jpoles], Abs[#1] > Abs[#2] &];
  nVoigt = Length[Jpoles];

  (* Residues of J at each Voigt pole *)
  JdenRPrime = D[JdenRPoly, ss];
  Jresidues = Table[
    Re[(JnumPoly /. ss -> Jpoles[[j]]) /
      (Jpoles[[j]] (JdenRPrime /. ss -> Jpoles[[j]]))],
    {j, nVoigt}];

  (* Extract Voigt parameters: alpha_i = spring, eta_i = dashpot (paper Eq. 13) *)
  etaVc = 1/Jresidues;
  alphaVc = -Jpoles/Jresidues;
  alphaVPa = alphaVc Fconv;
  etaVPa = etaVc Fconv;
  tauV = -1/Jpoles;

  (* Return all results *)
  <|
    "kE" -> kEval,
    "kf" -> kfval,
    "k2Re" -> k2Re,
    "k2Im" -> k2Im,
    "gamma" -> gammaPa,
    "alpha" -> alphaPa,
    "eta" -> eta0Pa,
    "prefactor" -> prefactor,
    "Fconv" -> Fconv,
    "nPoles" -> nPoles,
    "poles" -> poles9,
    "residues" -> residj,
    "nVoigt" -> nVoigt,
    "tauV" -> tauV,
    "alphaV" -> alphaVPa,
    "etaV" -> etaVPa
  |>
];

End[];

EndPackage[];
