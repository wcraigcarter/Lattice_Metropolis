(* Regenerates demo.png, the README header figure: a kinetic-tool anneal
   with attractive like-atom bonds, coarsening domains, and exchange with
   the gas at the free boundary.  Positive chemical potentials keep the
   vacancy content realistic (light decoration of the domain boundaries
   rather than wide vacancy channels). *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "LatticeGasMC.wl"}]];
SeedRandom[4];

demoConfiguration = latticeGasRandomConfiguration[120, 160, 8500, 8500];
{demoConfiguration, statistics} = latticeGasMetropolisSweep[
   demoConfiguration, 5*10^7,
   "energyNearest" -> -1.0, "energyNextNearest" -> -0.4,
   "inverseTemperature" -> 2.0, "muOrange" -> 0.8, "muBlue" -> 0.8,
   "exchangeProbability" -> 0.1];
Export[FileNameJoin[{DirectoryName[$InputFileName], "demo.png"}],
  latticeGasConfigurationPlot[demoConfiguration, ImageSize -> 480]];
Print["exported demo.png;  counts: ",
  latticeGasSpeciesCounts[demoConfiguration]];
