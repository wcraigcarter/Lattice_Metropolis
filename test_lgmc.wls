(* Validation + benchmarks for LatticeGasMC` *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "LatticeGasMC.wl"}]];

SeedRandom[7];

(* ---------- Test 1a: conservation & frame invariance, exchangeProbability = 0 *)
configuration = latticeGasRandomConfiguration[30, 40, 350, 350];
initialCounts = latticeGasSpeciesCounts[configuration];
{configuration2, statistics} = latticeGasMetropolisSweep[configuration, 2*10^6,
   "energyNearest" -> 0.6, "energyNextNearest" -> 0.25,
   "inverseTemperature" -> 1.5];
Print["T1a species conserved (no exchanges):  ",
  latticeGasSpeciesCounts[configuration2] === initialCounts];
Print["T1a frame still all zero:              ",
  Total[Abs[configuration2[[1]]]] + Total[Abs[configuration2[[-1]]]] +
    Total[Abs[configuration2[[All, 1]]]] +
    Total[Abs[configuration2[[All, -1]]]] === 0];
Print["T1a vacancy-switch acceptance:         ",
  N[statistics["acceptedVacancySwitches"]/
    statistics["attemptedVacancySwitches"]]];

(* ---------- Test 1b: single long kernel call, unequal frequencies,
   exchanges on: atom-count bookkeeping must match the move statistics ----- *)
configuration = latticeGasRandomConfiguration[200, 200, 12000, 12000];
initialCounts = latticeGasSpeciesCounts[configuration];
{configuration2, statistics} = latticeGasMetropolisSweep[configuration, 10^7,
   "energyNearest" -> 0.6, "energyNextNearest" -> 0.25,
   "inverseTemperature" -> 1.5, "muOrange" -> 0.3, "muBlue" -> -0.2,
   "exchangeProbability" -> 0.2,
   "vacancySwitchFrequency" -> 1.0, "occupiedSwitchFrequency" -> 0.2];
finalCounts = latticeGasSpeciesCounts[configuration2];
Print["T1b atom-count change == insertions - removals:  ",
  (finalCounts["orange"] + finalCounts["blue"]) -
    (initialCounts["orange"] + initialCounts["blue"]) ===
   statistics["acceptedInsertions"] - statistics["acceptedRemovals"]];
Print["T1b frame zero, values in {-1,0,1}:              ",
  Total[Abs[configuration2[[1]]]] + Total[Abs[configuration2[[-1]]]] +
    Total[Abs[configuration2[[All, 1]]]] +
    Total[Abs[configuration2[[All, -1]]]] === 0 &&
   SubsetQ[{-1, 0, 1}, Union[Flatten[configuration2]]]];
Print["T1b statistics: ", statistics];

(* ---------- Test 2: exact distribution on 2 x 3 lattice (729 states),
   UNEQUAL switch frequencies (equilibrium must be unchanged) -------------- *)
energyNearest = 0.6; energyNextNearest = 0.25; muOrange = 0.4; muBlue = -0.2;
inverseTemperature = 1.2; exchangeProbability = 0.5;
vacancySwitchFrequency = 1.0; occupiedSwitchFrequency = 0.35;
states = ArrayReshape[#, {2, 3}] & /@ Tuples[{-1, 0, 1}, 6];
paddedStates = ArrayPad[#, 1] & /@ states;
energies = latticeGasTotalEnergy[#, energyNearest, energyNextNearest] & /@
   paddedStates;
orangeCounts = Count[#, 1, {2}] & /@ states;
blueCounts = Count[#, -1, {2}] & /@ states;
weights = Exp[-inverseTemperature (energies - muOrange orangeCounts -
      muBlue blueCounts)];
exactProbabilities = weights/Total[weights];

nSamples = 100000; stepsBetweenSamples = 60;
configuration = paddedStates[[RandomInteger[{1, 729}]]];
sampleCounts = ConstantArray[0, 729];
Do[
  configuration = First@latticeGasMetropolisSweepCore[configuration,
     energyNearest, energyNextNearest, muOrange, muBlue,
     inverseTemperature, exchangeProbability,
     vacancySwitchFrequency, occupiedSwitchFrequency,
     stepsBetweenSamples];
  sampleCounts[[FromDigits[
      Flatten[configuration[[2 ;; -2, 2 ;; -2]]] + 1, 3] + 1]]++,
  {nSamples}];
sampledProbabilities = sampleCounts/N[nSamples];
totalVariation = Total[Abs[sampledProbabilities - exactProbabilities]]/2;
Print["T2 total-variation dist (expect ~0.03 noise floor, fail >> 0.1):  ",
  N[totalVariation]];
Print["T2 <E>   exact vs sampled:  ", N[exactProbabilities . energies],
  "  ", N[sampledProbabilities . energies]];
Print["T2 <nO>  exact vs sampled:  ", N[exactProbabilities . orangeCounts],
  "  ", N[sampledProbabilities . orangeCounts]];
Print["T2 <nB>  exact vs sampled:  ", N[exactProbabilities . blueCounts],
  "  ", N[sampledProbabilities . blueCounts]];

(* ---------- Test 3: 3 x 3 lattice moments, PURELY VACANCY-MEDIATED
   switches (occupiedSwitchFrequency = 0) + exchanges --------------------- *)
energyNearest = 0.5; energyNextNearest = 0.2; muOrange = 0.35; muBlue = -0.15;
inverseTemperature = 1.0; exchangeProbability = 0.5;
vacancySwitchFrequency = 1.0; occupiedSwitchFrequency = 0.0;
states = ArrayReshape[#, {3, 3}] & /@ Tuples[{-1, 0, 1}, 9];
paddedStates = ArrayPad[#, 1] & /@ states;
energies = latticeGasTotalEnergy[#, energyNearest, energyNextNearest] & /@
   paddedStates;
orangeCounts = Count[#, 1, {2}] & /@ states;
blueCounts = Count[#, -1, {2}] & /@ states;
weights = Exp[-inverseTemperature (energies - muOrange orangeCounts -
      muBlue blueCounts)];
exactProbabilities = weights/Total[weights];

nSamples = 100000; stepsBetweenSamples = 80;
configuration = paddedStates[[RandomInteger[{1, 3^9}]]];
energySum = 0.; orangeSum = 0; blueSum = 0;
Do[
  configuration = First@latticeGasMetropolisSweepCore[configuration,
     energyNearest, energyNextNearest, muOrange, muBlue,
     inverseTemperature, exchangeProbability,
     vacancySwitchFrequency, occupiedSwitchFrequency,
     stepsBetweenSamples];
  energySum += latticeGasTotalEnergy[configuration, energyNearest,
     energyNextNearest];
  orangeSum += Count[configuration, 1, {2}];
  blueSum += Count[configuration, -1, {2}],
  {nSamples}];
Print["T3 <E>   exact vs sampled:  ", N[exactProbabilities . energies],
  "  ", energySum/nSamples];
Print["T3 <nO>  exact vs sampled:  ", N[exactProbabilities . orangeCounts],
  "  ", orangeSum/N[nSamples]];
Print["T3 <nB>  exact vs sampled:  ", N[exactProbabilities . blueCounts],
  "  ", blueSum/N[nSamples]];

(* ---------- Test 4: benchmarks at 500 x 500 ---------- *)
bigConfiguration = latticeGasRandomConfiguration[500, 500, 100000, 100000];
timePerAttempt = First@AbsoluteTiming[
     latticeGasMetropolisSweepCore[bigConfiguration,
      0.6, 0.25, 0., 0., 1.5, 0., 1., 1., 10^7]]/10^7;
Print["T4 equal frequencies:               ",
  NumberForm[10^9 timePerAttempt, 3], " ns/attempt"];
(* dilute vacancies (0.5%), vacancy-mediated only: every attempt lands on a
   vacancy thanks to the list; uniform site selection would waste 99.5% *)
diluteConfiguration = latticeGasRandomConfiguration[500, 500, 124375, 124375];
timePerAttempt = First@AbsoluteTiming[
     latticeGasMetropolisSweepCore[diluteConfiguration,
      0.6, 0.25, 0., 0., 1.5, 0., 1., 0., 10^7]]/10^7;
Print["T4 0.5% vacancies, vacancy-only:    ",
  NumberForm[10^9 timePerAttempt, 3], " ns/attempt"];
timePerAttempt = First@AbsoluteTiming[
     latticeGasMetropolisSweepCore[bigConfiguration,
      0.6, 0.25, 0.3, -0.2, 1.5, 0.1, 1., 1., 10^7]]/10^7;
Print["T4 with 10% exchanges:              ",
  NumberForm[10^9 timePerAttempt, 3], " ns/attempt"];

(* ---------- Test 6: BULK-EXCHANGE tool, exact distribution on 2 x 3 ------
   must reproduce the same semi-grand ensemble as the surface-exchange tool *)
energyNearest = 0.6; energyNextNearest = 0.25; muOrange = 0.4; muBlue = -0.2;
inverseTemperature = 1.2;
states = ArrayReshape[#, {2, 3}] & /@ Tuples[{-1, 0, 1}, 6];
paddedStates = ArrayPad[#, 1] & /@ states;
energies = latticeGasTotalEnergy[#, energyNearest, energyNextNearest] & /@
   paddedStates;
orangeCounts = Count[#, 1, {2}] & /@ states;
blueCounts = Count[#, -1, {2}] & /@ states;
weights = Exp[-inverseTemperature (energies - muOrange orangeCounts -
      muBlue blueCounts)];
exactProbabilities = weights/Total[weights];

nSamples = 100000; stepsBetweenSamples = 40;
configuration = paddedStates[[RandomInteger[{1, 729}]]];
sampleCounts = ConstantArray[0, 729];
Do[
  configuration = First@latticeGasBulkExchangeSweepCore[configuration,
     energyNearest, energyNextNearest, muOrange, muBlue,
     inverseTemperature, stepsBetweenSamples];
  sampleCounts[[FromDigits[
      Flatten[configuration[[2 ;; -2, 2 ;; -2]]] + 1, 3] + 1]]++,
  {nSamples}];
sampledProbabilities = sampleCounts/N[nSamples];
totalVariation = Total[Abs[sampledProbabilities - exactProbabilities]]/2;
Print["T6 bulk tool, total-variation dist (expect ~0.03 noise floor):  ",
  N[totalVariation]];
Print["T6 <E>   exact vs sampled:  ", N[exactProbabilities . energies],
  "  ", N[sampledProbabilities . energies]];
Print["T6 <nO>  exact vs sampled:  ", N[exactProbabilities . orangeCounts],
  "  ", N[sampledProbabilities . orangeCounts]];
Print["T6 <nB>  exact vs sampled:  ", N[exactProbabilities . blueCounts],
  "  ", N[sampledProbabilities . blueCounts]];

(* ---------- Test 7: bulk tool on 3 x 3, moments ---------- *)
energyNearest = 0.5; energyNextNearest = 0.2; muOrange = 0.35; muBlue = -0.15;
inverseTemperature = 1.0;
states = ArrayReshape[#, {3, 3}] & /@ Tuples[{-1, 0, 1}, 9];
paddedStates = ArrayPad[#, 1] & /@ states;
energies = latticeGasTotalEnergy[#, energyNearest, energyNextNearest] & /@
   paddedStates;
orangeCounts = Count[#, 1, {2}] & /@ states;
blueCounts = Count[#, -1, {2}] & /@ states;
weights = Exp[-inverseTemperature (energies - muOrange orangeCounts -
      muBlue blueCounts)];
exactProbabilities = weights/Total[weights];

nSamples = 100000; stepsBetweenSamples = 60;
configuration = paddedStates[[RandomInteger[{1, 3^9}]]];
energySum = 0.; orangeSum = 0; blueSum = 0;
Do[
  configuration = First@latticeGasBulkExchangeSweepCore[configuration,
     energyNearest, energyNextNearest, muOrange, muBlue,
     inverseTemperature, stepsBetweenSamples];
  energySum += latticeGasTotalEnergy[configuration, energyNearest,
     energyNextNearest];
  orangeSum += Count[configuration, 1, {2}];
  blueSum += Count[configuration, -1, {2}],
  {nSamples}];
Print["T7 <E>   exact vs sampled:  ", N[exactProbabilities . energies],
  "  ", energySum/nSamples];
Print["T7 <nO>  exact vs sampled:  ", N[exactProbabilities . orangeCounts],
  "  ", orangeSum/N[nSamples]];
Print["T7 <nB>  exact vs sampled:  ", N[exactProbabilities . blueCounts],
  "  ", blueSum/N[nSamples]];

(* ---------- Test 8: bulk tool benchmark at 500 x 500 ---------- *)
bigConfiguration = latticeGasRandomConfiguration[500, 500, 100000, 100000];
timePerAttempt = First@AbsoluteTiming[
     latticeGasBulkExchangeSweepCore[bigConfiguration,
      0.6, 0.25, 0.3, -0.2, 1.5, 10^7]]/10^7;
Print["T8 bulk-exchange tool:              ",
  NumberForm[10^9 timePerAttempt, 3], " ns/attempt"];

(* ---------- Test 9: order-parameter maps on perfect synthetic patterns ----
   each perfect phase must give |own map| ~ 1 (interior) and others ~ 0 --- *)
smoothingRadius = 4; margin = 2 smoothingRadius + 2;
interior[map_] := map[[margin ;; -margin, margin ;; -margin]];
mapNames = {"composition", "checkerboard", "stripeVertical",
   "stripeHorizontal"};
patterns = <|
   "composition" -> Table[1, {i, 40}, {j, 40}],
   "checkerboard" -> Table[(-1)^(i + j), {i, 40}, {j, 40}],
   "stripeVertical" -> Table[(-1)^j, {i, 40}, {j, 40}],
   "stripeHorizontal" -> Table[(-1)^i, {i, 40}, {j, 40}]|>;
allPassed = True;
Do[
  maps = latticeGasOrderParameterMaps[ArrayPad[patterns[patternName], 1],
    smoothingRadius];
  ownSignal = Min[Abs[interior[maps[patternName]]]];
  crossTalk = Max[Table[Max[Abs[interior[maps[otherName]]]],
     {otherName, DeleteCases[mapNames, patternName]}]];
  If[ownSignal < 0.999 || crossTalk > 0.02,
   allPassed = False;
   Print["T9 FAIL for ", patternName, ": own ", ownSignal,
    ", crosstalk ", crossTalk]],
  {patternName, mapNames}];
maps = latticeGasOrderParameterMaps[
   ArrayPad[Table[(-1)^(i + j), {i, 40}, {j, 40}], 1], smoothingRadius];
Print["T9 perfect patterns: own map ~ +-1, others ~ 0:   ", allPassed];
Print["T9 occupancy of dense pattern ~ 1:                ",
  Min[interior[maps["occupancy"]]] > 0.999];

(* ---------- Test 10: phase map classifies synthetic quadrants ---------- *)
(* 80 x 80: four quadrants = solid orange | checkerboard |
                             vertical stripes | horizontal stripes *)
quadrants = ArrayPad[
   Table[Which[
     i <= 40 && j <= 40, 1,
     i <= 40, (-1)^(i + j),
     j <= 40, (-1)^j,
     True, (-1)^i],
    {i, 80}, {j, 80}], 1];
phaseMap = latticeGasPhaseMap[quadrants, 4];
quadrantCenter[rowRange_, columnRange_] :=
  Union[Flatten[phaseMap[[rowRange, columnRange]]]];
Print["T10 quadrant centers classified correctly:  ",
  quadrantCenter[12 ;; 28, 12 ;; 28] === {1} &&
   MatchQ[quadrantCenter[12 ;; 28, 52 ;; 68], {3} | {4}] &&
   MatchQ[quadrantCenter[52 ;; 68, 12 ;; 28], {5} | {6}] &&
   MatchQ[quadrantCenter[52 ;; 68, 52 ;; 68], {7} | {8}]];
(* random configuration is mostly disordered; empty region is vacancy *)
disorderMap = latticeGasPhaseMap[
   ArrayPad[RandomChoice[{1, -1}, {60, 60}], 1], 4];
emptyMap = latticeGasPhaseMap[ArrayPad[Table[0, {40}, {40}], 1], 4];
Print["T10 random -> mostly disordered (code 9):    ",
  Count[Flatten[disorderMap[[10 ;; -10, 10 ;; -10]]], 9] >
   0.9 Length[Flatten[disorderMap[[10 ;; -10, 10 ;; -10]]]]];
Print["T10 empty lattice -> vacancy (code 0):       ",
  Union[Flatten[emptyMap]] === {0}];

(* ---------- Test 5: plotting ---------- *)
demoConfiguration = latticeGasRandomConfiguration[120, 160, 5000, 5000];
{demoConfiguration, statistics} = latticeGasMetropolisSweep[
   demoConfiguration, 5*10^7,
   "energyNearest" -> -1.0, "energyNextNearest" -> -0.4,
   "inverseTemperature" -> 2.0, "muOrange" -> -0.5, "muBlue" -> -0.5,
   "exchangeProbability" -> 0.1];
Export[FileNameJoin[{DirectoryName[$InputFileName], "demo.png"}],
  latticeGasConfigurationPlot[demoConfiguration, ImageSize -> 480]];
Print["T5 demo plot exported;  counts after anneal: ",
  latticeGasSpeciesCounts[demoConfiguration]];
