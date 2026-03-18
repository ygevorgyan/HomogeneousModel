(* ::Package:: *)

(* ================================================================= *)
(* VisualizeLayers.wl -- 2D cross-section of a multilayered body     *)
(* ================================================================= *)
(* Reads ModelInput.wl, draws a to-scale sector cross-section with   *)
(* numbered layers on the left, property table on the right.          *)
(* Works for any number of layers.                                    *)
(* ================================================================= *)

(* ----------------------------------------------------------------- *)
(* Section 0: Load model                                              *)
(* ----------------------------------------------------------------- *)
If[StringQ[$InputFileName] && $InputFileName =!= "",
  SetDirectory[DirectoryName[$InputFileName]],
  SetDirectory[NotebookDirectory[]]
];
Get["ModelInput.wl"];

rk   = $ModelConfig["rk"];
muk  = $ModelConfig["muk"];
nLayers = Length[rk];

(* ----------------------------------------------------------------- *)
(* Section 1: Derived quantities                                      *)
(* ----------------------------------------------------------------- *)

rkKm = N[rk / 1000];
rSurface = Last[rkKm];

(* Inner boundary radii: 0 for layer 1, then previous outer radii *)
rInner = Prepend[Most[rkKm], 0.];
thicknessKm = rkKm - rInner;

(* Liquid detection *)
liquidQ = Table[N[muk[[i]]] < 1, {i, nLayers}];

(* Layer names: use config if provided, otherwise generic "Layer i (solid/liquid)" *)
layerNames = If[KeyExistsQ[$ModelConfig, "layerNames"],
  $ModelConfig["layerNames"],
  Table[
    "Layer " <> ToString[i] <> If[liquidQ[[i]], " (liquid)", " (solid)"],
    {i, nLayers}
  ]
];

(* Colors: LightTemperatureMap palette, warm at center -> cool at surface *)
(* Map to [0.2, 0.8] to avoid extreme saturation at both ends *)
layerColors = Table[
  ColorData["LightTemperatureMap"][0.8 - 0.6 * (i - 1) / (nLayers - 1)],
  {i, nLayers}
];

(* Normalized radii (0 to 1) -- true to scale *)
normOuter = rkKm / rSurface;
normInner = rInner / rSurface;

(* Thin-layer exaggeration: enforce minimum visual fraction per layer *)
minVisualFrac = 0.03;
thicknessExaggerated = False;
Module[{fracs, deficit, thinIdx, thickIdx, totalSteal, stealFrac, cumulative},
  fracs = normOuter - normInner;
  thinIdx  = Flatten[Position[fracs, _?(# < minVisualFrac &)]];
  thickIdx = Complement[Range[nLayers], thinIdx];
  If[Length[thinIdx] > 0 && Length[thickIdx] > 0,
    deficit = Total[(minVisualFrac - fracs[[#]]) & /@ thinIdx];
    totalSteal = Total[fracs[[thickIdx]]];
    (* Steal proportionally from thick layers *)
    stealFrac = deficit / totalSteal;
    Do[fracs[[j]] = fracs[[j]] (1 - stealFrac), {j, thickIdx}];
    Do[fracs[[j]] = minVisualFrac, {j, thinIdx}];
    (* Rebuild normOuter/normInner from adjusted fractions *)
    cumulative = Accumulate[fracs];
    normOuter = cumulative;
    normInner = Prepend[Most[cumulative], 0.];
    thicknessExaggerated = True;
  ];
];

(* ----------------------------------------------------------------- *)
(* Section 2: Build cross-section graphic                             *)
(* ----------------------------------------------------------------- *)

(* Sector angular range: 60 deg symmetric around vertical (Pi/2) *)
thetaMin = Pi / 3;      (* 60 deg *)
thetaMax = 2 Pi / 3;    (* 120 deg *)

(* Draw layers as filled sectors, outermost first *)
sectorLayers = Table[
  With[{i = nLayers + 1 - k},
    {EdgeForm[{Thickness[0.003], GrayLevel[0.3]}],
     layerColors[[i]],
     Disk[{0, 0}, normOuter[[i]], {thetaMin, thetaMax}]}
  ],
  {k, 1, nLayers}
];

(* Boundary arcs *)
boundaryArcs = Table[
  {Dashed, GrayLevel[0.4], Thickness[0.003],
   Circle[{0, 0}, normOuter[[i]], {thetaMin, thetaMax}]},
  {i, 1, nLayers - 1}
];

(* Radial edge lines *)
edgeLines = {Thickness[0.003], GrayLevel[0.3],
  Line[{{0, 0}, {Cos[thetaMin], Sin[thetaMin]}}],
  Line[{{0, 0}, {Cos[thetaMax], Sin[thetaMax]}}]};

(* Common style *)
lblSize = 22;
lblFont = FontFamily -> "Times";

(* Number formatting: one decimal place, no unit suffix *)
fmtKm[x_] := ToString[NumberForm[N[Round[x, 0.1]], {Infinity, 1}]];

(* Layer number labels at midpoint of each layer's left boundary *)
layerNumberLabels = Table[
  Module[{rMid, xMid, yMid},
    rMid = (normInner[[i]] + normOuter[[i]]) / 2;
    (* Midpoint on the left radial edge at this radius *)
    xMid = rMid * Cos[thetaMax];
    yMid = rMid * Sin[thetaMax];
    Text[Style[ToString[i], lblSize - 2, GrayLevel[0.15], lblFont],
      {xMid - 0.06, yMid}, {1, 0}]   (* offset left so it sits outside the edge *)
  ],
  {i, 1, nLayers}
];

(* Build the cross-section graphic *)
sectorGraphic = Graphics[
  {sectorLayers, boundaryArcs, edgeLines, layerNumberLabels},
  PlotRange -> All,
  PlotRangePadding -> {{Scaled[0.08], Scaled[0.02]}, {Scaled[0.05], Scaled[0.02]}},
  AspectRatio -> Automatic,
  ImageSize -> 300
];

(* Build the layer info listing: "i --- Name  thickness km (radius km)" *)
tblSize = lblSize - 4;
layerListing = Column[
  Table[
    Style[
      ToString[i] <> " \[LongDash] " <> layerNames[[i]] <> "   " <>
        fmtKm[thicknessKm[[i]]] <> " km (" <> fmtKm[rkKm[[i]]] <> " km)",
      tblSize, FontFamily -> "Times"],
    {i, 1, nLayers}
  ],
  Spacings -> 1
];

(* Title: note exaggeration if applied *)
titleSuffix = If[thicknessExaggerated,
  " (visual thickness exaggerated)", " (to scale)"];
titleText = Style[
  "Layered body \[LongDash] cross-section" <> titleSuffix,
  Bold, lblSize, FontFamily -> "Times"
];

(* Combine: title above, cross-section + table side-by-side *)
fig = Column[{
  titleText,
  Grid[{{sectorGraphic, layerListing}},
    Spacings -> {4, 0},
    Alignment -> {Center, Center}]
}, Alignment -> Center, Spacings -> 2];

(* ----------------------------------------------------------------- *)
(* Section 3: Export                                                   *)
(* ----------------------------------------------------------------- *)
outputPath = FileNameJoin[{Directory[], "output", "layer_cross_section.pdf"}];
Export[outputPath, fig, ImageSize -> 1200];
Print["Cross-section exported to: ", outputPath];
