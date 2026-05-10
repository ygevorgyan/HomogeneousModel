(* ::Package:: *)

(* ================================================================= *)
(* ModelInputStratLiquid.wl -- stratified liquid core demonstration   *)
(* ================================================================= *)
(* Five-layer Maxwell body in which the liquid Fe core is supplied   *)
(* as three stratified liquid sub-layers, with two solid silicate   *)
(* mantle layers above:                                                *)
(*                                                                    *)
(*   Layer 1: Liquid inner Fe core   (radii 0--2000 km, rho 12000)    *)
(*   Layer 2: Liquid outer Fe core   (radii 2000--3500 km, rho 10000) *)
(*   Layer 3: Liquid Fe outer shell  (radii 3500--3700 km, rho 8000)  *)
(*   Layer 4: Lower silicate mantle  (radii 3700--5500 km, rho 4500)  *)
(*   Layer 5: Upper silicate mantle  (radii 5500--6371 km, rho 3300)  *)
(*                                                                    *)
(* The mergeAdjacentLiquidLayers pass should fuse layers 1-3 into a   *)
(* single liquid layer with outer radius 3700 km and a volume-weighted *)
(* mean density of 10008.77 kg/m^3 (preserves total mass exactly).    *)
(* ================================================================= *)

$ModelConfig = <|
  "rk"     -> {2000000, 3500000, 3700000, 5500000, 6371000},
  "rhok"   -> {12000, 10000, 8000, 4500, 3300},
  "muk"    -> {10^-6, 10^-6, 10^-6, 80 10^9, 65 10^9},
  "etak"   -> {1, 1, 1, 10^21, 10^21},
  "omega"  -> 2 Pi / (1 * 86400),
  "layerNames" -> {"Liquid inner Fe", "Liquid outer Fe",
                    "Liquid Fe shell", "Lower mantle",
                    "Upper mantle"}
|>;

$AnalysisConfig = <|
  "workingPrecision"  -> 50,
  "gridLogMin"        -> -17,
  "gridLogMax"        -> -5,
  "gridPoints"        -> 500,
  "omegaInterestMin"  -> 10^-9,
  "omegaInterestMax"  -> 10^-3,
  "dominantThreshold" -> 1/50
|>;
