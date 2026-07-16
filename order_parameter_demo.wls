(* Demonstrate order-parameter maps on three simulated microstructures *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "LatticeGasMC.wl"}]];
SeedRandom[11];

(* 1. decomposition: like atoms attract *)
decomposition = latticeGasRandomConfiguration[200, 200, 16000, 16000];
{decomposition, s} = latticeGasBulkExchangeSweep[decomposition, 3*10^7,
   "energyNearest" -> -1.0, "energyNextNearest" -> -0.4,
   "inverseTemperature" -> 1.5, "muOrange" -> 0., "muBlue" -> 0.];

(* 2. checkerboard: unlike nearest neighbors preferred *)
checkerboard = latticeGasRandomConfiguration[200, 200, 16000, 16000];
{checkerboard, s} = latticeGasBulkExchangeSweep[checkerboard, 3*10^7,
   "energyNearest" -> 1.0, "energyNextNearest" -> -0.2,
   "inverseTemperature" -> 1.5, "muOrange" -> 1.0, "muBlue" -> 1.0];

(* 3. stripes: unlike nearest AND unlike next-nearest preferred,
      energyNextNearest > energyNearest/2 *)
stripes = latticeGasRandomConfiguration[200, 200, 16000, 16000];
{stripes, s} = latticeGasBulkExchangeSweep[stripes, 3*10^7,
   "energyNearest" -> 1.0, "energyNextNearest" -> 1.0,
   "inverseTemperature" -> 1.5, "muOrange" -> 1.0, "muBlue" -> 1.0];

figure = Column[{
    Labeled[latticeGasOrderParameterPlot[decomposition, 8],
     Style["decomposition  (energyNearest = -1, energyNextNearest = -0.4)",
      14], Top],
    Labeled[latticeGasOrderParameterPlot[checkerboard, 8],
     Style["checkerboard  (energyNearest = +1, energyNextNearest = -0.2)",
      14], Top],
    Labeled[latticeGasOrderParameterPlot[stripes, 8],
     Style["stripes  (energyNearest = +1, energyNextNearest = +1)", 14],
     Top]},
   Spacings -> 2];
Export[FileNameJoin[{DirectoryName[$InputFileName],
   "order_parameters.png"}], figure, ImageResolution -> 120];
Print["exported order_parameters.png"];

(* one-panel-per-microstructure phase maps sharing a single legend *)
phaseMapFigure = Legended[
   Grid[{
     {Labeled[latticeGasConfigurationPlot[decomposition,
        ImageSize -> 260], "decomposition", Top],
      Labeled[latticeGasConfigurationPlot[checkerboard,
        ImageSize -> 260], "checkerboard", Top],
      Labeled[latticeGasConfigurationPlot[stripes,
        ImageSize -> 260], "stripes", Top]},
     {latticeGasPhaseMapPlot[decomposition, 8, 0.4, "Legend" -> False,
       ImageSize -> 260],
      latticeGasPhaseMapPlot[checkerboard, 8, 0.4, "Legend" -> False,
       ImageSize -> 260],
      latticeGasPhaseMapPlot[stripes, 8, 0.4, "Legend" -> False,
       ImageSize -> 260]}},
    Spacings -> {1, 1}],
   latticeGasPhaseMapLegend[]];
Export[FileNameJoin[{DirectoryName[$InputFileName],
   "phase_maps.png"}], phaseMapFigure, ImageResolution -> 120];
Print["exported phase_maps.png"];
