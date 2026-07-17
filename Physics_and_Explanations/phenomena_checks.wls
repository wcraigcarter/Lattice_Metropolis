(* Verifies the quantitative claims in Phenomena.md and regenerates the
   figures in figures/.  Run:  wolframscript -file phenomena_checks.wls  *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "LatticeGasMC.wl"}]];
SeedRandom[17];
figuresDirectory = FileNameJoin[{DirectoryName[$InputFileName], "figures"}];
If[! DirectoryQ[figuresDirectory], CreateDirectory[figuresDirectory]];

checkerboardOrderParameter[configuration_] :=
  Module[{realSites = configuration[[2 ;; -2, 2 ;; -2]], nRows, nColumns},
   {nRows, nColumns} = Dimensions[realSites];
   Total[Table[(-1)^(i + j), {i, nRows}, {j, nColumns}] realSites, 2]/
    (nRows nColumns)];

(* ---- 1. checkerboard order-disorder transition: |OP| vs temperature ---- *)
Print["--- 1. order-disorder transition (energyNearest = +1, energyNextNearest = -0.2) ---"];
configuration = latticeGasRandomConfiguration[64, 64, 2048, 2048];
orderParameterData = Table[
   {configuration, statistics} =
    latticeGasBulkExchangeSweep[configuration, 10^7,
     "energyNearest" -> 1.0, "energyNextNearest" -> -0.2,
     "muOrange" -> 1.5, "muBlue" -> 1.5,
     "inverseTemperature" -> 1/temperature];
   {temperature, Abs[checkerboardOrderParameter[configuration]]},
   {temperature, Range[3.4, 0.6, -0.2]}];
Print[TableForm[orderParameterData,
   TableHeadings -> {None, {"T", "|checkerboard OP|"}}]];
Export[FileNameJoin[{figuresDirectory, "order_parameter_vs_temperature.png"}],
  ListLinePlot[orderParameterData, PlotMarkers -> Automatic,
   Frame -> True,
   FrameLabel -> {"temperature  (units of |energyNearest|, kB = 1)",
     "|checkerboard order parameter|"},
   PlotRange -> {{0.5, 3.45}, {-0.02, 1.05}},
   GridLines -> {{2.269}, None},
   PlotLabel -> "gridline: NN-only Ising Tc = 2.269", ImageSize -> 520]];

(* ---- 2. curvature-driven shrinkage of a circular antiphase domain ------ *)
Print["--- 2. antiphase disk shrinkage (Allen-Cahn: area ~ linear in time) ---"];
nSide = 128; diskRadius = 36;
checkerboardPattern = Table[(-1)^(i + j), {i, nSide}, {j, nSide}];
configuration = ArrayPad[
   Table[checkerboardPattern[[i, j]] If[
      (i - (nSide + 1)/2.)^2 + (j - (nSide + 1)/2.)^2 <= diskRadius^2,
      -1, 1], {i, nSide}, {j, nSide}], 1];
minorityRegistryArea[config_] :=
  Count[Flatten[checkerboardPattern config[[2 ;; -2, 2 ;; -2]]], -1];
attemptsPerStep = 10^6; nSites = nSide^2;
shrinkageData = {{0., minorityRegistryArea[configuration]}};
Do[
  {configuration, statistics} =
   latticeGasBulkExchangeSweep[configuration, attemptsPerStep,
    "energyNearest" -> 1.0, "energyNextNearest" -> -0.2,
    "muOrange" -> 1.5, "muBlue" -> 1.5, "inverseTemperature" -> 1.5];
  AppendTo[shrinkageData,
   {step attemptsPerStep/N[nSites], minorityRegistryArea[configuration]}];
  If[shrinkageData[[-1, 2]] == 0, Break[]],
  {step, 1, 200}];
linearFit = Fit[shrinkageData[[;; -2]], {1, t}, t];
Print["initial area (pi r^2 = ", Round[Pi diskRadius^2], "): ",
  shrinkageData[[1, 2]], ";  gone after ", shrinkageData[[-1, 1]],
  " sweeps;  linear fit: ", linearFit];
Export[FileNameJoin[{figuresDirectory, "antiphase_disk_shrinkage.png"}],
  Show[
   ListPlot[shrinkageData, Frame -> True,
    FrameLabel -> {"time  (sweeps = attempts per site)",
      "antiphase-domain area  (sites)"},
    PlotRange -> All, ImageSize -> 520],
   Plot[linearFit, {t, 0, shrinkageData[[-1, 1]]},
    PlotStyle -> Directive[Dashed, Gray]]]];

(* ---- 3. crystal-vapor: droplet growth / coexistence / sublimation ------ *)
Print["--- 3. single-species droplet vs muOrange (T = 0 coexistence at -2.8) ---"];
dropletConfiguration[radius_] := ArrayPad[
   Table[If[(i - 50.5)^2 + (j - 50.5)^2 <= radius^2, 1, 0],
    {i, 100}, {j, 100}], 1];
dropletCurves = Table[
   configuration = dropletConfiguration[25];
   curve = Table[
     {configuration, statistics} =
      latticeGasBulkExchangeSweep[configuration, 5*10^5,
       "energyNearest" -> -1.0, "energyNextNearest" -> -0.4,
       "muOrange" -> muOrange, "muBlue" -> -6.0,
       "inverseTemperature" -> 1.5];
     {sample 0.5, latticeGasSpeciesCounts[configuration]["orange"]},
     {sample, 1, 40}];
   Print["  muOrange = ", muOrange, ":  atoms ",
    Round[Pi 25^2], " -> ", curve[[-1, 2]]];
   curve,
   {muOrange, {-2.4, -2.8, -3.2}}];
Export[FileNameJoin[{figuresDirectory, "droplet_growth_vs_mu.png"}],
  ListLinePlot[dropletCurves, Frame -> True,
   FrameLabel -> {"attempts  (millions)", "atoms in the system"},
   PlotLegends -> {"muOrange = -2.4  (supersaturated: growth)",
     "muOrange = -2.8  (T = 0 coexistence)",
     "muOrange = -3.2  (undersaturated: sublimation)"},
   PlotRange -> All, ImageSize -> 520]];

(* ---- 4. equilibrium vacancy concentration vs mu (Arrhenius) ------------ *)
Print["--- 4. bulk vacancy concentration vs muOrange (predict exp[-beta(mu + 5.6)]) ---"];
Do[
  configuration = latticeGasRandomConfiguration[100, 100, 10^4, 0];
  {configuration, statistics} =
   latticeGasBulkExchangeSweep[configuration, 10^7,
    "energyNearest" -> -1.0, "energyNextNearest" -> -0.4,
    "muOrange" -> muOrange, "muBlue" -> -6.0, "inverseTemperature" -> 1.0];
  vacancySamples = Table[
    {configuration, statistics} =
     latticeGasBulkExchangeSweep[configuration, 2*10^5,
      "energyNearest" -> -1.0, "energyNextNearest" -> -0.4,
      "muOrange" -> muOrange, "muBlue" -> -6.0,
      "inverseTemperature" -> 1.0];
    latticeGasSpeciesCounts[configuration]["vacancy"], {20}];
  Print["  muOrange = ", muOrange, ":  measured <nVacancy> = ",
   N[Mean[vacancySamples]], ",  predicted 10^4 exp[-(mu+5.6)] = ",
   NumberForm[10^4 Exp[-(muOrange + 5.6)], 4]],
  {muOrange, {1.0, 1.5, 2.0, 2.5}}];

(* ---- 5. first-order transition and hysteresis in delta mu -------------- *)
Print["--- 5. composition hysteresis vs deltaMu = muOrange - muBlue ---"];
deltaMuValues = Range[-2.0, 2.0, 0.2];
orangeFraction[config_] := Module[{counts = latticeGasSpeciesCounts[config]},
   counts["orange"]/N[counts["orange"] + counts["blue"]]];
sweepBranch[deltas_] := Module[{config, branch = {}},
   config = latticeGasRandomConfiguration[64, 64, 1850, 1850];
   Do[
    {config, statistics} = latticeGasBulkExchangeSweep[config, 2*10^6,
      "energyNearest" -> -1.0, "energyNextNearest" -> -0.4,
      "muOrange" -> 1. + deltaMu/2., "muBlue" -> 1. - deltaMu/2.,
      "inverseTemperature" -> 1.0];
    AppendTo[branch, {deltaMu, orangeFraction[config]}],
    {deltaMu, deltas}];
   branch];
upBranch = sweepBranch[deltaMuValues];
downBranch = sweepBranch[Reverse[deltaMuValues]];
Print["  up-sweep   fractions: ", NumberForm[#, 2] & /@ upBranch[[All, 2]]];
Print["  down-sweep fractions: ", NumberForm[#, 2] & /@ downBranch[[All, 2]]];
Export[FileNameJoin[{figuresDirectory, "composition_hysteresis.png"}],
  ListLinePlot[{upBranch, downBranch}, PlotMarkers -> Automatic,
   Frame -> True,
   FrameLabel -> {"deltaMu = muOrange - muBlue", "orange fraction of atoms"},
   PlotLegends -> {"sweep up", "sweep down"},
   PlotRange -> {{-2.1, 2.1}, {-0.05, 1.05}}, ImageSize -> 520]];

(* ---- 6. ergodicity: vacancy-free lattice with direct exchange off ------ *)
Print["--- 6. frozen dynamics without vacancies ---"];
denseConfiguration = ArrayPad[RandomChoice[{1, -1}, {30, 30}], 1];
{frozenConfiguration, frozenStatistics} =
  latticeGasMetropolisSweep[denseConfiguration, 10^6,
   "energyNearest" -> -1.0, "occupiedSwitchFrequency" -> 0.];
Print["  0 vacancies, occupiedSwitchFrequency = 0:  configuration unchanged = ",
  frozenConfiguration === denseConfiguration, ";  statistics = ",
  frozenStatistics];
denseConfiguration[[15, 15]] = 0;
{thawedConfiguration, thawedStatistics} =
  latticeGasMetropolisSweep[denseConfiguration, 10^6,
   "energyNearest" -> -1.0, "occupiedSwitchFrequency" -> 0.];
Print["  1 vacancy:  acceptedVacancySwitches = ",
  thawedStatistics["acceptedVacancySwitches"]];
