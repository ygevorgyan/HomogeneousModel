(* ::Package:: *)
(* ================================================================= *)
(* RunHomogeneous.wl — Standalone Equivalent Homogeneous Body        *)
(* ================================================================= *)
(* Computes the equivalent homogeneous Voigt rheology from a layered *)
(* model defined in ModelInput.wl, without running the full          *)
(* sensitivity pipeline.                                              *)
(*                                                                     *)
(* Run with: wolframscript -file RunHomogeneous.wl                   *)
(* ================================================================= *)

t0global = AbsoluteTime[];


(* ================================================================= *)
(* Section 0: Configuration & Module Loading                          *)
(* ================================================================= *)

If[StringQ[$InputFileName] && $InputFileName =!= "",
  SetDirectory[DirectoryName[$InputFileName]],
  SetDirectory[NotebookDirectory[]]
];

Get["src/LayeredBodyConstants.wl"];
Get["src/YMatrix.wl"];
Get["src/NLayerJacobian.wl"];
Get["src/HomogeneousRheology.wl"];

(* Load user-editable configuration *)
Get["ModelInput.wl"];

(* Unpack analysis settings *)
wp = $AnalysisConfig["workingPrecision"];
$NLwp = wp;
gridLogMin = $AnalysisConfig["gridLogMin"];
gridLogMax = $AnalysisConfig["gridLogMax"];
gridPoints = $AnalysisConfig["gridPoints"];
omegaIntMin = $AnalysisConfig["omegaInterestMin"];
omegaIntMax = $AnalysisConfig["omegaInterestMax"];
dominantThreshold = $AnalysisConfig["dominantThreshold"];

(* Unpack physical model *)
rk0 = SetPrecision[$ModelConfig["rk"], wp];
rhok0 = SetPrecision[$ModelConfig["rhok"], wp];
muk0 = SetPrecision[$ModelConfig["muk"], wp];
etak0 = SetPrecision[$ModelConfig["etak"], wp];
omega = SetPrecision[$ModelConfig["omega"], wp];

(* Auto-detect layer structure *)
nLayers = Length[rk0];

If[Length[rhok0] != nLayers || Length[muk0] != nLayers || Length[etak0] != nLayers,
  Print["ERROR: All parameter arrays in $ModelConfig must have the same length."];
  Quit[1]];
If[nLayers < 2,
  Print["ERROR: Need at least 2 layers."]; Quit[1]];

Print["\n========================================"];
Print["  ", ToString[nLayers], "-Layer Body -> Homogeneous Rheology"];
Print["========================================\n"];

Print["Configuration:"];
Print["  Working precision = ", wp];
Print["  omega = ", ScientificForm[N[omega], 6], " rad/s"];
Print["  Layers: ", nLayers];
Print["  Grid: log10|s| in [", gridLogMin, ", ", gridLogMax, "], ", gridPoints, " points"];
Print["  Interest range: omega in [", ScientificForm[N[omegaIntMin], 3],
  ", ", ScientificForm[N[omegaIntMax], 3], "] rad/s"];
Print["  Dominant mode threshold: ", NumberForm[100 N[dominantThreshold], 3], "%"];
Print["  Modules loaded: LayeredBodyConstants, YMatrix, NLayerJacobian, HomogeneousRheology"];


(* ================================================================= *)
(* Section 1: Compute Homogeneous Rheology                            *)
(* ================================================================= *)

Print["\n--- Computing equivalent homogeneous body ---"];
t0 = AbsoluteTime[];
result = computeHomogeneous[rk0, rhok0, muk0, etak0, omega,
  "WorkingPrecision" -> wp,
  "GridLogMin" -> gridLogMin,
  "GridLogMax" -> gridLogMax,
  "GridPoints" -> gridPoints];
t1 = AbsoluteTime[];
Print["Computed in ", NumberForm[t1 - t0, 4], " seconds."];


(* ================================================================= *)
(* Section 2: Print Results                                            *)
(* ================================================================= *)

Print["\n========================================================="];
Print["  EQUIVALENT HOMOGENEOUS BODY"];
Print["  Rheology: Generalized Voigt (Gevorgyan et al. 2023)"];
Print["========================================================="];

Print["\n--- Love numbers ---"];
Print["  kE       = ", ScientificForm[result["kE"], 6]];
Print["  kf       = ", ScientificForm[result["kf"], 6]];
Print["  Re[k2]   = ", ScientificForm[result["k2Re"], 6]];
Print["  Im[k2]   = ", ScientificForm[result["k2Im"], 6]];

Print["\n--- Rheological parameters ---"];
Print["  gamma    = ", fmtSci[result["gamma"]/10^9, 4], " GPa"];
Print["  alpha    = ", fmtSci[result["alpha"]/10^9, 4], " GPa"];
Print["  eta      = ", fmtSci[result["eta"], 4], " Pa s"];
Print["  nPoles   = ", result["nPoles"]];
Print["  nVoigt   = ", result["nVoigt"]];

yrSec = SetPrecision[365.25 24 3600, wp];
nP = result["nPoles"];

Print["\n--- Secular relaxation modes (sorted by |s|, fastest first) ---"];
Print[StringPadRight["  j", 6],
  StringPadRight["s_j (1/s)", 16],
  StringPadRight["tau_j (s)", 16],
  StringPadRight["tau_j (yr)", 14],
  StringPadRight["residue r_j", 16]];
Print["  ", StringJoin[Table["-", 66]]];

(* poles are already sorted by decreasing |s| from computeHomogeneous *)
Do[
  sj = result["poles"][[j]];
  rj = result["residues"][[j]];
  tauj = -1/sj;
  Print[StringPadRight["  " <> ToString[j], 6],
    StringPadRight[fmtSci[sj, 4], 16],
    StringPadRight[fmtSci[tauj, 4], 16],
    StringPadRight[fmtSci[tauj/yrSec, 4], 14],
    StringPadRight[fmtSci[rj, 4], 16]],
  {j, nP}];

Print["\n--- Voigt elements (sorted by tau, slowest first) ---"];
Print[StringPadRight["  i", 6],
  StringPadRight["tau (s)", 16],
  StringPadRight["tau (yr)", 14],
  StringPadRight["alpha (GPa)", 16],
  StringPadRight["eta (Pa s)", 16]];
Print["  ", StringJoin[Table["-", 66]]];

tauSorted = Sort[result["tauV"], #1 > #2 &];
sortIdx = Ordering[result["tauV"], All, #1 > #2 &];

Do[
  jj = sortIdx[[i]];
  Print[StringPadRight["  " <> ToString[i], 6],
    StringPadRight[fmtSci[tauSorted[[i]], 4], 16],
    StringPadRight[fmtSci[tauSorted[[i]]/yrSec, 4], 14],
    StringPadRight[fmtSci[result["alphaV"][[jj]]/10^9, 4], 16],
    StringPadRight[fmtSci[result["etaV"][[jj]], 4], 16]],
  {i, Length[tauSorted]}];


(* ================================================================= *)
(* Section 2b: Dominant Mode Identification                           *)
(* ================================================================= *)

(* Compliance-space quantities (used here and in Section 3 plots) *)
prefactor = result["prefactor"];
Fconv = result["Fconv"];
gammac = prefactor / result["kf"];
alphac = prefactor / result["kE"] - gammac;
eta0c = result["eta"] / Fconv;
nV = result["nVoigt"];
muVc = result["alphaV"] / Fconv;
etaVc = result["etaV"] / Fconv;

(* Sweep the user-specified frequency interest range to find each
   Voigt element's peak fractional contribution to -Im[J].
   A mode is dominant if its max fractional contribution anywhere
   in the interest range exceeds the threshold. *)
nSweep = 200;
omegaSweep = Exp[Subdivide[Log[N[omegaIntMin]], Log[N[omegaIntMax]], nSweep]];

(* For each sweep frequency, compute per-element and total -Im[J] *)
(* (excluding the spring 1/alpha which is real and doesn't contribute) *)
maxFracEl = Table[0., nV];
Do[
  sv = I omS;
  imDash = -Im[1/(eta0c sv)];
  imEls = Table[-Im[1/(muVc[[j]] + etaVc[[j]] sv)], {j, nV}];
  imTot = imDash + Total[imEls];
  If[imTot > 0,
    Do[
      frac = imEls[[j]] / imTot;
      If[frac > maxFracEl[[j]], maxFracEl[[j]] = frac],
      {j, nV}]],
  {omS, omegaSweep}];

(* Also compute contribution at the tidal frequency for display *)
sOmega = I omega;
imJelTidal = Table[-Im[1/(muVc[[j]] + etaVc[[j]] sOmega)], {j, nV}];
imJdashTidal = -Im[1/(eta0c sOmega)];
imJtotalTidal = imJdashTidal + Total[imJelTidal];
fracElTidal = imJelTidal / imJtotalTidal;
fracDashTidal = imJdashTidal / imJtotalTidal;

(* Sort by max broadband contribution (descending) *)
contribOrder = Ordering[maxFracEl, All, #1 > #2 &];

Print["\n--- Mode contributions to dissipation ---"];
Print["  Interest range: [", ScientificForm[N[omegaIntMin], 3],
  ", ", ScientificForm[N[omegaIntMax], 3], "] rad/s"];
Print[StringPadRight["  Mode", 8],
  StringPadRight["tau (yr)", 14],
  StringPadRight["@ tidal", 10],
  StringPadRight["Max in range", 14],
  "Dominant?"];
Print["  ", StringJoin[Table["-", 60]]];
Do[
  jj = contribOrder[[i]];
  isDom = If[maxFracEl[[jj]] > dominantThreshold, " *", ""];
  Print[StringPadRight["  V" <> ToString[jj], 8],
    StringPadRight[fmtSci[result["tauV"][[jj]]/yrSec, 4], 14],
    StringPadRight[ToString[NumberForm[100 fracElTidal[[jj]], {4, 1}]] <> "%", 10],
    StringPadRight[ToString[NumberForm[100 maxFracEl[[jj]], {4, 1}]] <> "%", 14],
    isDom],
  {i, nV}];
Print[StringPadRight["  dashpot", 8],
  StringPadRight["Inf", 14],
  ToString[NumberForm[100 fracDashTidal, {4, 1}]] <> "%",
  "  (always kept)"];

(* Select dominant modes: max fractional contribution in interest range > threshold *)
dominantIdx = Select[Range[nV], maxFracEl[[#]] > dominantThreshold &];
nDominant = Length[dominantIdx];

Print["\nDominant modes (max contribution > ",
  ToString[NumberForm[100 N[dominantThreshold], 3]],
  "% in interest range): ", nDominant, " of ", nV];
Print["  Indices: ", dominantIdx];
Print["  Combined tidal-freq contribution: ",
  NumberForm[100 Total[fracElTidal[[dominantIdx]]], {4, 1}], "%"];

(* Build reduced Voigt model *)
JReduced[sv_] := 1/alphac + 1/(eta0c sv) +
  Sum[1/(muVc[[dominantIdx[[j]]]] + etaVc[[dominantIdx[[j]]]] sv),
    {j, nDominant}];
k2Reduced[sv_] := prefactor JReduced[sv] / (1 + gammac JReduced[sv]);

(* Evaluate reduced model at tidal frequency *)
k2RedVal = k2Reduced[sOmega];
Print["\n--- Reduced model (", nDominant, " Voigt elements) ---"];
Print["  Re[k2]   = ", ScientificForm[Re[k2RedVal], 6],
  "  (full: ", ScientificForm[result["k2Re"], 6], ")"];
Print["  Im[k2]   = ", ScientificForm[Im[k2RedVal], 6],
  "  (full: ", ScientificForm[result["k2Im"], 6], ")"];
Print["  Rel. error Re[k2]: ",
  NumberForm[100 Abs[(Re[k2RedVal] - result["k2Re"])/result["k2Re"]], {3, 2}], "%"];
Print["  Rel. error Im[k2]: ",
  NumberForm[100 Abs[(Im[k2RedVal] - result["k2Im"])/result["k2Im"]], {3, 2}], "%"];

(* Print reduced Voigt element table *)
Print["\n--- Reduced Voigt elements ---"];
Print[StringPadRight["  i", 6],
  StringPadRight["tau (s)", 16],
  StringPadRight["tau (yr)", 14],
  StringPadRight["alpha (GPa)", 16],
  StringPadRight["eta (Pa s)", 16]];
Print["  ", StringJoin[Table["-", 66]]];

domTauSorted = Sort[result["tauV"][[dominantIdx]], #1 > #2 &];
domSortIdx = dominantIdx[[Ordering[result["tauV"][[dominantIdx]], All, #1 > #2 &]]];

Do[
  jj = domSortIdx[[i]];
  Print[StringPadRight["  " <> ToString[i], 6],
    StringPadRight[fmtSci[domTauSorted[[i]], 4], 16],
    StringPadRight[fmtSci[domTauSorted[[i]]/yrSec, 4], 14],
    StringPadRight[fmtSci[result["alphaV"][[jj]]/10^9, 4], 16],
    StringPadRight[fmtSci[result["etaV"][[jj]], 4], 16]],
  {i, Length[domTauSorted]}];


(* ================================================================= *)
(* Section 3: Love Number Visualization                               *)
(* ================================================================= *)

Print["\n--- Generating Love number plots ---"];

(* Frequency grid matching the interest range, 200 points log-spaced *)
omegaPlot = Exp[Subdivide[Log[N[omegaIntMin]], Log[N[omegaIntMax]], 200]];

(* Evaluate k2 at each frequency via the N-layer numerical formula *)
k2Multi = Table[
  nLayerK2Num[SetPrecision[I omegaP, wp], rk0, rhok0, muk0, etak0],
  {omegaP, omegaPlot}];

(* Reconstruct k2 from the extracted Voigt parameters *)
(* Compliance-space vars (prefactor, gammac, etc.) defined in Section 2b *)
JVoigt[sv_] := 1/alphac + 1/(eta0c sv) +
  Sum[1/(muVc[[j]] + etaVc[[j]] sv), {j, nV}];
k2Recon[sv_] := prefactor JVoigt[sv] / (1 + gammac JVoigt[sv]);

k2VRecon = Table[k2Recon[I omegaP], {omegaP, omegaPlot}];

(* Reduced model reconstruction (dominant modes only) *)
k2VReduced = Table[k2Reduced[I omegaP], {omegaP, omegaPlot}];

(* Plot 1: -Im[k2] — full Voigt (all modes) *)
plot1 = ListLogLogPlot[{
  Transpose[{omegaPlot, -Im[k2Multi]}],
  Transpose[{omegaPlot, -Im[k2VRecon]}]},
  PlotStyle -> {
    {Blue, Thick},
    {Red, Dashed, Thick}},
  Joined -> True,
  PlotRange -> All,
  AxesLabel -> {"\[Omega] (rad/s)", "-Im[k2]"},
  PlotLabel -> "Dissipation (all " <> ToString[nV] <> " Voigt elements)",
  PlotLegends -> {"N-layer", "Full Voigt"},
  BaseStyle -> {FontFamily -> "Times", FontSize -> 16},
  LabelStyle -> Directive[FontSize -> 18],
  TicksStyle -> Directive[FontSize -> 14],
  ImageSize -> Large];
Print[plot1]

(* Plot 2: -Im[k2] — reduced Voigt (dominant modes) *)
plot2 = ListLogLogPlot[{
  Transpose[{omegaPlot, -Im[k2Multi]}],
  Transpose[{omegaPlot, -Im[k2VReduced]}]},
  PlotStyle -> {
    {Blue, Thick},
    {Darker[Green], DotDashed, Thick}},
  Joined -> True,
  PlotRange -> All,
  AxesLabel -> {"\[Omega] (rad/s)", "-Im[k2]"},
  PlotLabel -> "Dissipation (dominant " <> ToString[nDominant] <> " Voigt elements)",
  PlotLegends -> {"N-layer", "Reduced Voigt"},
  BaseStyle -> {FontFamily -> "Times", FontSize -> 16},
  LabelStyle -> Directive[FontSize -> 18],
  TicksStyle -> Directive[FontSize -> 14],
  ImageSize -> Large];
Print[plot2]

(* Plot 3: Re[k2] — full Voigt (all modes) *)
plot3 = ListLogLinearPlot[{
  Transpose[{omegaPlot, Re[k2Multi]}],
  Transpose[{omegaPlot, Re[k2VRecon]}]},
  PlotStyle -> {
    {Blue, Thick},
    {Red, Dashed, Thick}},
  Joined -> True,
  PlotRange -> All,
  AxesLabel -> {"\[Omega] (rad/s)", "Re[k2]"},
  PlotLabel -> "Elastic response (all " <> ToString[nV] <> " Voigt elements)",
  PlotLegends -> {"N-layer", "Full Voigt"},
  BaseStyle -> {FontFamily -> "Times", FontSize -> 16},
  LabelStyle -> Directive[FontSize -> 18],
  TicksStyle -> Directive[FontSize -> 14],
  ImageSize -> Large];
Print[plot3]

(* Plot 4: Re[k2] — reduced Voigt (dominant modes) *)
plot4 = ListLogLinearPlot[{
  Transpose[{omegaPlot, Re[k2Multi]}],
  Transpose[{omegaPlot, Re[k2VReduced]}]},
  PlotStyle -> {
    {Blue, Thick},
    {Darker[Green], DotDashed, Thick}},
  Joined -> True,
  PlotRange -> All,
  AxesLabel -> {"\[Omega] (rad/s)", "Re[k2]"},
  PlotLabel -> "Elastic response (dominant " <> ToString[nDominant] <> " Voigt elements)",
  PlotLegends -> {"N-layer", "Reduced Voigt"},
  BaseStyle -> {FontFamily -> "Times", FontSize -> 16},
  LabelStyle -> Directive[FontSize -> 18],
  TicksStyle -> Directive[FontSize -> 14],
  ImageSize -> Large];
Print[plot4]

(* Export plots to PDF *)
outDir = FileNameJoin[{Directory[], "output"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir]];

Export[FileNameJoin[{outDir, "imag_k2_full.pdf"}], plot1];
Export[FileNameJoin[{outDir, "imag_k2_reduced.pdf"}], plot2];
Export[FileNameJoin[{outDir, "real_k2_full.pdf"}], plot3];
Export[FileNameJoin[{outDir, "real_k2_reduced.pdf"}], plot4];

Print["\nPlots exported to output/:"];
Print["  imag_k2_full.pdf, imag_k2_reduced.pdf"];
Print["  real_k2_full.pdf, real_k2_reduced.pdf"];


(* ================================================================= *)
(* Section 4: Export Model Parameters                                  *)
(* ================================================================= *)

Print["\n--- Exporting model parameters ---"];

(* Build output association *)
outputParams = <|
  "loveNumbers" -> <|
    "kE"   -> N[result["kE"]],
    "kf"   -> N[result["kf"]],
    "k2Re" -> N[result["k2Re"]],
    "k2Im" -> N[result["k2Im"]]
  |>,
  "rheology" -> <|
    "gamma_Pa" -> N[result["gamma"]],
    "alpha_Pa" -> N[result["alpha"]],
    "eta_Pas"  -> N[result["eta"]],
    "nPoles"   -> result["nPoles"],
    "nVoigt"   -> result["nVoigt"]
  |>,
  "poles" -> Table[<|
    "s_1ps"     -> N[result["poles"][[j]]],
    "tau_s"     -> N[-1/result["poles"][[j]]],
    "tau_yr"    -> N[-1/result["poles"][[j]] / yrSec],
    "residue"   -> N[result["residues"][[j]]]
  |>, {j, result["nPoles"]}],
  "voigtElements" -> Table[<|
    "alpha_Pa" -> N[result["alphaV"][[j]]],
    "eta_Pas"  -> N[result["etaV"][[j]]],
    "tau_s"    -> N[result["tauV"][[j]]],
    "tau_yr"   -> N[result["tauV"][[j]] / yrSec]
  |>, {j, result["nVoigt"]}],
  "dominantModes" -> <|
    "threshold"  -> N[dominantThreshold],
    "indices"    -> dominantIdx,
    "nDominant"  -> nDominant,
    "k2Re_reduced" -> N[Re[k2RedVal]],
    "k2Im_reduced" -> N[Im[k2RedVal]]
  |>,
  "inputConfig" -> <|
    "nLayers" -> nLayers,
    "omega_radps" -> N[omega],
    "workingPrecision" -> wp
  |>
|>;

(* Export JSON *)
Export[FileNameJoin[{outDir, "homogeneous_params.json"}], outputParams, "JSON"];

(* Export Wolfram .wl *)
Export[FileNameJoin[{outDir, "homogeneous_params.wl"}], outputParams];

Print["  homogeneous_params.json"];
Print["  homogeneous_params.wl"];

tfinal = AbsoluteTime[];
Print["\nTotal runtime: ", NumberForm[tfinal - t0global, 4], " seconds."];

Print["\n========================================"];
Print["  Done."];
Print["========================================\n"];
