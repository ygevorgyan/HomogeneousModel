(* ::Package:: *)

(* ================================================================= *)
(* ModelInputEnceladus.wl -- 3-layer Enceladus model                 *)
(* ================================================================= *)
(* Based on Cadek et al. (2016), Thomas et al. (2016).               *)
(* Tests: 3-layer model, existing layer-name auto-detection for      *)
(* nLayers==3.                                                        *)
(*                                                                    *)
(*   Layer 1: Silicate core    (solid)                                *)
(*   Layer 2: Subsurface ocean (liquid)                               *)
(*   Layer 3: Ice shell        (solid)                                *)
(* ================================================================= *)

$ModelConfig = <|
  "rk"   -> {192000, 230000, 252100},                     (* m *)
  "rhok" -> {2400, 1000, 930},                             (* kg/m^3 *)
  "muk"  -> {50 10^9, 1/10, 33 10^8},                     (* Pa; layer 2 liquid *)
  "etak" -> {10^19, 1, 10^14},                             (* Pa s *)
  "omega" -> 2 Pi / (13702/10000 * 86400),                 (* rad/s; 1.3702-day period *)
  "layerNames" -> {"Silicate core", "Subsurface ocean", "Ice shell"}
|>;

$AnalysisConfig = <|
  "workingPrecision"  -> 50,
  "gridLogMin"        -> -18,
  "gridLogMax"        -> -3,
  "gridPoints"        -> 2000,
  "omegaInterestMin"  -> 10^-10,
  "omegaInterestMax"  -> 10^-3,
  "dominantThreshold" -> 1/100
|>;
