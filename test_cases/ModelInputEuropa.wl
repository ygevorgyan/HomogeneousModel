(* ::Package:: *)

(* ================================================================= *)
(* ModelInputEuropa.wl -- 4-layer Europa model                       *)
(* ================================================================= *)
(* Based on Wahr et al. (2006), Hussmann & Spohn (2004).             *)
(* Tests: 4-layer model with two liquid layers (core + ocean).       *)
(*                                                                    *)
(*   Layer 1: Iron core       (liquid)                                *)
(*   Layer 2: Silicate mantle (solid)                                 *)
(*   Layer 3: Subsurface ocean (liquid)                               *)
(*   Layer 4: Ice shell        (solid)                                *)
(* ================================================================= *)

$ModelConfig = <|
  "rk"   -> {600000, 1400000, 1540000, 1560800},         (* m *)
  "rhok" -> {5150, 3300, 1000, 930},                      (* kg/m^3 *)
  "muk"  -> {1/10, 65 10^9, 1/10, 35 10^8},              (* Pa; layers 1,3 liquid *)
  "etak" -> {1, 10^20, 1, 10^14},                         (* Pa s *)
  "omega" -> 2 Pi / (35512/10000 * 86400),                (* rad/s; 3.5512-day period *)
  "layerNames" -> {"Iron core", "Silicate mantle", "Subsurface ocean", "Ice shell"}
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
