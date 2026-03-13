(* ::Package:: *)

(* ================================================================= *)
(* ModelInput.wl -- 5-layer Moon model                                *)
(* ================================================================= *)
(* Matsuyama et al. (2016) / Briaud (2023).                           *)
(*                                                                    *)
(*   Layer 1: Inner core       (solid iron)                           *)
(*   Layer 2: Outer core       (liquid iron)                          *)
(*   Layer 3: Lower mantle     (partial melt zone)                    *)
(*   Layer 4: Upper mantle     (solid silicate)                       *)
(*   Layer 5: Crust            (solid silicate)                       *)
(* ================================================================= *)

$ModelConfig = <|
  "rk"   -> {213517, 325000, 500000, 1697150, 1737150},    (* m *)
  "rhok" -> {772016/100, 6700, 3800, 3356, 2735},          (* kg/m^3 *)
  "muk"  -> {40 10^9, 1/10, 60 10^9, 625 10^8, 15 10^9},  (* Pa; layer 2 liquid *)
  "etak" -> {10^21, 1, 5 10^16, 10^21, 10^23},             (* Pa s *)
  "omega" -> 2 Pi / (2732/100 * 86400),                    (* rad/s; 27.32-day period *)
  "layerNames" -> {"Inner core", "Outer core", "Lower mantle", "Upper mantle", "Crust"}
|>;

$AnalysisConfig = <|
  "workingPrecision"  -> 50,
  "gridLogMin"        -> -20,
  "gridLogMax"        -> -3,
  "gridPoints"        -> 2000,
  "omegaInterestMin"  -> 10^-10,
  "omegaInterestMax"  -> 10^-5,
  "dominantThreshold" -> 1/50
|>;
