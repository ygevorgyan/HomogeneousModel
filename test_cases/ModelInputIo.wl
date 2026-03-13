(* ::Package:: *)

(* ================================================================= *)
(* ModelInputIo.wl -- 4-layer Io model                               *)
(* ================================================================= *)
(* Based on Segatz et al. (1988), Bierson & Nimmo (2016).            *)
(* Tests: 4-layer model with one liquid layer, contrasting            *)
(* viscosities in solid layers.                                       *)
(*                                                                    *)
(*   Layer 1: Fe-FeS core     (liquid)                                *)
(*   Layer 2: Mantle           (solid)                                *)
(*   Layer 3: Asthenosphere    (solid, low viscosity)                 *)
(*   Layer 4: Lithosphere      (solid, high viscosity)                *)
(* ================================================================= *)

$ModelConfig = <|
  "rk"   -> {952000, 1591000, 1791600, 1821600},         (* m *)
  "rhok" -> {5150, 3244, 3244, 3244},                     (* kg/m^3 *)
  "muk"  -> {1/10, 60 10^9, 60 10^9, 65 10^9},           (* Pa; layer 1 liquid *)
  "etak" -> {1, 10^20, 10^16, 10^23},                     (* Pa s *)
  "omega" -> 2 Pi / (17691/10000 * 86400),                 (* rad/s; 1.7691-day period *)
  "layerNames" -> {"Fe-FeS core", "Mantle", "Asthenosphere", "Lithosphere"}
|>;

$AnalysisConfig = <|
  "workingPrecision"  -> 50,
  "gridLogMin"        -> -20,
  "gridLogMax"        -> -3,
  "gridPoints"        -> 2000,
  "omegaInterestMin"  -> 10^-10,
  "omegaInterestMax"  -> 10^-3,
  "dominantThreshold" -> 1/100
|>;
