Homogeneous equivalent body — incompressible model
===================================================

Computes Love numbers and the equivalent homogeneous Voigt rheology
for a multilayered incompressible body under tidal forcing.


Prerequisites
-------------
Wolfram Engine (free) or Mathematica, version 12.0+.
Verify that wolframscript is available:

    wolframscript -version

If it is not on your PATH but Mathematica is installed, call it
directly, e.g.:

    /Applications/Mathematica.app/Contents/MacOS/wolframscript -version


Configuration
-------------
Edit ModelInput.wl to define the body.  It contains two blocks:

  $ModelConfig      radii, densities, shear moduli, viscosities,
                    tidal frequency, optional layer names
  $AnalysisConfig   precision, pole-search grid, frequency range,
                    dominant-mode threshold

All layer arrays (rk, rhok, muk, etak) must have the same length
and are ordered from center (layer 1) to surface (layer N).
Liquid layers are indicated by setting mu < 1 Pa.

An optional "layerNames" key can provide custom labels for the
visualization (e.g., {"Inner core", "Outer core", ...}).  If
omitted, layers are labelled generically as "Layer i (solid/liquid)".

Pre-configured test cases for the Moon, Europa, Io, and Enceladus
are in test_cases/.  To use one, copy it over ModelInput.wl:

    cp test_cases/ModelInputEuropa.wl ModelInput.wl


Running
-------
All commands are run from the homogeneous_model/ directory.

1. Equivalent homogeneous body computation:

    wolframscript -file RunHomogeneous.wl

   Outputs Love numbers (kE, kf, k2), rheological parameters
   (gamma, alpha, eta), secular relaxation modes, Voigt elements,
   dominant mode analysis, and four comparison plots.
   Runtime: under 1 minute.

2. Layer cross-section figure:

    wolframscript -file VisualizeLayers.wl

   Produces output/layer_cross_section.pdf.  Runs in a few seconds.


Output
------
Terminal:  tables of Love numbers, rheological parameters, secular
           modes, Voigt elements, dominant mode analysis.

Files in output/:
  homogeneous_params.wl     Voigt parameters (Wolfram format)
  homogeneous_params.json   Voigt parameters (JSON)
  real_k2_full.pdf          Re[k2]: N-layer vs full Voigt
  real_k2_reduced.pdf       Re[k2]: N-layer vs reduced Voigt
  imag_k2_full.pdf          -Im[k2]: N-layer vs full Voigt
  imag_k2_reduced.pdf       -Im[k2]: N-layer vs reduced Voigt
  layer_cross_section.pdf   Layer cross-section diagram


File layout
-----------
ModelInput.wl         User configuration
RunHomogeneous.wl     Main driver script
VisualizeLayers.wl    Layer cross-section generator
src/                  Core modules:
  LayeredBodyConstants.wl   Physical constants, gravity, formatting
  YMatrix.wl                6x6 fundamental matrix and propagator
  NLayerJacobian.wl         N-layer propagator, k2(s)
  HomogeneousRheology.wl    Equivalent homogeneous body computation
test_cases/           Alternative model configurations:
  ModelInputMoon.wl         5-layer Moon
  ModelInputEuropa.wl       4-layer Europa
  ModelInputIo.wl           4-layer Io
  ModelInputEnceladus.wl    3-layer Enceladus


Troubleshooting
---------------
- FindRoot::brdig warnings are harmless.
- If poles are missed (fewer modes than expected), widen
  gridLogMin/gridLogMax or increase gridPoints in $AnalysisConfig.
- If plots show separation between N-layer and full Voigt curves,
  increase workingPrecision to 80 or 100.


References
----------
Gevorgyan, Matsuyama & Ragazzo (2023), MNRAS 523, 1822.

Sabadini, Vermeersen & Cambiotti (2016), Global Dynamics of the
Earth: Applications of Viscoelastic Relaxation Theory to Solid-Earth
and Planetary Geophysics, 2nd ed., Springer.
