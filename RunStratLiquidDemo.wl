(* ::Package:: *)
(* Demonstration that adjacent liquid sub-layers in a Maxwell body are
   merged automatically into one equivalent liquid layer, and that the
   resulting k_2 matches the result of supplying the merged layer
   directly. *)

If[StringQ[$InputFileName] && $InputFileName =!= "",
  SetDirectory[DirectoryName[$InputFileName]],
  SetDirectory[NotebookDirectory[]]];

Get["src/LayeredBodyConstants.wl"];
Get["src/YMatrix.wl"];
Get["src/NLayerJacobian.wl"];
Get["src/HomogeneousRheology.wl"];

Get["test_cases/ModelInputStratLiquid.wl"];

wp = $AnalysisConfig["workingPrecision"];
$NLwp = wp;

rk    = SetPrecision[$ModelConfig["rk"],   wp];
rhok  = SetPrecision[$ModelConfig["rhok"], wp];
muk   = SetPrecision[$ModelConfig["muk"],  wp];
etak  = SetPrecision[$ModelConfig["etak"], wp];
omega = SetPrecision[$ModelConfig["omega"], wp];

Print["\n=========================================================="];
Print["  Stratified-liquid demonstration (Maxwell pipeline)"];
Print["=========================================================="];
Print["  Input layers (", Length[rk], "):"];
Do[
  Print["    ", j, ". ",
    StringPadRight[$ModelConfig["layerNames"][[j]], 22],
    "  r = ", N[rk[[j]]/1000], " km    rho = ",
    N[rhok[[j]]], " kg/m^3"],
  {j, Length[rk]}];

(* --- Show what mergeAdjacentLiquidLayers does --- *)
Print["\n  Merging adjacent liquid sub-layers ..."];
{rkM, rhokM, mukM, etakM} = mergeAdjacentLiquidLayers[
  rk, rhok, muk, etak];
Print["  Merged layers (", Length[rkM], "):"];
Do[
  Print["    ", j, ". r = ", N[rkM[[j]]/1000], " km    rho = ",
    NumberForm[N[rhokM[[j]]], 6], " kg/m^3    mu = ",
    ScientificForm[N[mukM[[j]]], 4], " Pa"],
  {j, Length[rkM]}];

(* --- Run computeHomogeneous on the stratified input (auto-merge) --- *)
Print["\n  computeHomogeneous on the stratified input ",
  "(merge on by default) ..."];
t0 = AbsoluteTime[];
resultStrat = computeHomogeneous[rk, rhok, muk, etak, omega,
  "WorkingPrecision" -> wp,
  "GridLogMin"  -> $AnalysisConfig["gridLogMin"],
  "GridLogMax"  -> $AnalysisConfig["gridLogMax"],
  "GridPoints"  -> $AnalysisConfig["gridPoints"]];
t1 = AbsoluteTime[];
Print["  done in ", NumberForm[t1 - t0, 4], " s"];

(* --- Run computeHomogeneous on the pre-merged input --- *)
Print["\n  computeHomogeneous on the same body supplied ",
  "as a single liquid layer ..."];
t0 = AbsoluteTime[];
resultMerged = computeHomogeneous[rkM, rhokM, mukM, etakM, omega,
  "WorkingPrecision" -> wp,
  "GridLogMin"  -> $AnalysisConfig["gridLogMin"],
  "GridLogMax"  -> $AnalysisConfig["gridLogMax"],
  "GridPoints"  -> $AnalysisConfig["gridPoints"]];
t1 = AbsoluteTime[];
Print["  done in ", NumberForm[t1 - t0, 4], " s"];

(* --- Compare --- *)
Print["\n  Comparison (stratified-as-supplied vs.\\ pre-merged):"];
Print["    kE         strat = ", ScientificForm[N[resultStrat["kE"]], 8]];
Print["               merged= ", ScientificForm[N[resultMerged["kE"]], 8]];
Print["    kf         strat = ", ScientificForm[N[resultStrat["kf"]], 8]];
Print["               merged= ", ScientificForm[N[resultMerged["kf"]], 8]];
Print["    Re[k2]     strat = ",
  ScientificForm[N[resultStrat["k2Re"]], 8]];
Print["               merged= ",
  ScientificForm[N[resultMerged["k2Re"]], 8]];
Print["    Im[k2]     strat = ",
  ScientificForm[N[resultStrat["k2Im"]], 8]];
Print["               merged= ",
  ScientificForm[N[resultMerged["k2Im"]], 8]];
Print["    nPoles     strat = ", resultStrat["nPoles"]];
Print["               merged= ", resultMerged["nPoles"]];

Print["\n    relative differences:"];
Print["      kE     : ", ScientificForm[
  N[(resultStrat["kE"] - resultMerged["kE"])/resultMerged["kE"]], 4]];
Print["      kf     : ", ScientificForm[
  N[(resultStrat["kf"] - resultMerged["kf"])/resultMerged["kf"]], 4]];
Print["      Re[k2] : ", ScientificForm[
  N[(resultStrat["k2Re"] - resultMerged["k2Re"])/resultMerged["k2Re"]], 4]];
Print["      Im[k2] : ", ScientificForm[
  N[(resultStrat["k2Im"] - resultMerged["k2Im"])/resultMerged["k2Im"]], 4]];
