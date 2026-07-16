(* ::Package:: *)

(* LatticeGasMC: Metropolis Monte Carlo for a two-species + vacancy lattice gas
   on an N x M square lattice with free (rectangular) boundaries.

   Encoding:  +1 = orange atom, -1 = blue atom, 0 = vacancy.
   Pair energy of neighboring sites i,j:  e * s_i * s_j
     (+e if same species, -e if unlike, 0 if either is a vacancy),
   with e = energyNearest for the 4 first (N,S,E,W) neighbors and
        e = energyNextNearest for the 4 second (diagonal) neighbors.

   Configurations are stored PADDED: an (N+2) x (M+2) integer matrix whose
   outer frame is permanently 0.  The zero frame is the free surface: bonds
   into it contribute nothing, and no move ever writes to it.

   Moves (one attempt per step):
     - with probability exchangeProbability: surface exchange with an
       infinite gas reservoir.  Pick a site uniformly from the outermost ring
       of real sites, pick a species s in {+1,-1} uniformly;
       if the site holds s, propose removal to the gas   (dE = mu_s - s*h);
       if the site is vacant, propose insertion of s     (dE = s*h - mu_s);
       otherwise reject.  h is the site enthalpy (weighted 8-neighbor sum).
       This proposal scheme is symmetric, so Metropolis acceptance
       min(1, exp(-beta dE)) satisfies detailed balance with the reservoir
       (semi-grand ensemble, weight
        exp[-beta (E - muOrange nOrange - muBlue nBlue)]).
     - otherwise, a switch attempt in one of two disjoint classes, chosen
       with probability proportional to
         vacancySwitchFrequency  * (number of vacancies)   vs
         occupiedSwitchFrequency * (number of atoms),
       so that every vacancy and every atom carries its own attempt
       frequency regardless of composition (kinetic-Monte-Carlo-style):
         * vacancy switch: a vacancy (drawn uniformly from a maintained list
           of vacancy positions) and one of its 8 neighbors, proceeding only
           if the neighbor is an atom -- vacancy-mediated diffusion;
         * occupied switch: an atom (drawn uniformly from a maintained list
           of atom positions) and one of its 8 neighbors, proceeding only if
           the neighbor is an atom of the other species -- direct exchange.
       Both classes use Metropolis acceptance on
         dE = (sa - sb) ((hb - jab sa) - (ha - jab sb)),
       where jab is the coupling of the a-b bond.  Both classes leave the
       vacancy and atom counts unchanged, so the class-selection weights are
       identical before and after any switch and detailed balance is exact.
       Only the RATIO of the two frequencies matters.  Setting
       occupiedSwitchFrequency = 0 gives purely vacancy-mediated (physical)
       kinetics; equal frequencies are equivalent to picking any real site
       uniformly.  With exchangeProbability = 0 both species counts are
       conserved (and at least one vacancy is then needed for
       vacancy-mediated dynamics to be ergodic).

   The vacancy/atom position lists are rebuilt by one O(N M) scan at the
   start of each kernel call and maintained in O(1) per accepted move.

   TWO SIMULATION TOOLS share this configuration format and these utilities:

   1. latticeGasMetropolisSweep -- the KINETIC tool described above:
      switches (vacancy-mediated and/or direct) plus optional exchange with
      the gas restricted to the surface.  Use it when the transport pathway
      matters (coarsening, diffusion, surface growth).

   2. latticeGasBulkExchangeSweep -- the PHASE-DIAGRAM tool: every attempt
      picks a real site uniformly (anywhere, not just the surface) and
      proposes changing its identity to one of the other two states
      {orange, blue, vacancy}, chosen with probability 1/2 each (a
      symmetric proposal).  Metropolis acceptance on
        dE = (new - old) * h - (mu_new - mu_old),   mu(vacancy) = 0,
      which covers insertion, removal, and direct orange <-> blue
      transmutation in one formula.  Nothing is conserved; the stationary
      distribution is the same semi-grand ensemble
      exp[-beta (E - muOrange nOrange - muBlue nBlue)] as tool 1 with
      exchanges on.  Equilibration is dramatically faster because no mass
      transport is needed -- the right tool for mapping phase diagrams
      (decomposition, order/disorder).  The two tools may be alternated on
      the same configuration; each is exact on its own, so any mixture is.

   Requires M >= 3 (per-run statistics ride home in matrix rows). *)

BeginPackage["LatticeGasMC`"];

latticeGasRandomConfiguration::usage =
  "latticeGasRandomConfiguration[nRows, nColumns, nOrangeAtoms, nBlueAtoms] \
gives a zero-padded (nRows+2) x (nColumns+2) configuration with the \
requested numbers of +1's and -1's placed uniformly at random among the \
real sites.  Counts are clipped (orange first, then blue) so that at least \
one vacancy always remains.";

latticeGasMetropolisSweep::usage =
  "latticeGasMetropolisSweep[configuration, nAttempts, opts] performs \
nAttempts Metropolis move attempts on the padded configuration and returns \
{newConfiguration, statistics}.  Options: \"energyNearest\", \
\"energyNextNearest\", \"muOrange\", \"muBlue\", \"inverseTemperature\", \
\"exchangeProbability\", \"vacancySwitchFrequency\", \
\"occupiedSwitchFrequency\".";

latticeGasMetropolisSweepCore::usage =
  "latticeGasMetropolisSweepCore[configuration, energyNearest, \
energyNextNearest, muOrange, muBlue, inverseTemperature, \
exchangeProbability, vacancySwitchFrequency, occupiedSwitchFrequency, \
nAttempts] is the low-overhead fixed-argument form of \
latticeGasMetropolisSweep; returns {newConfiguration, statistics}.";

latticeGasBulkExchangeSweep::usage =
  "latticeGasBulkExchangeSweep[configuration, nAttempts, opts] performs \
nAttempts bulk identity-change attempts (any real site may become orange, \
blue, or vacant, exchanging atoms with the reservoir regardless of \
position) and returns {newConfiguration, statistics}.  Samples the same \
semi-grand ensemble as latticeGasMetropolisSweep with exchanges on, but \
equilibrates much faster; intended for phase-diagram work where kinetics \
do not matter.  Options: \"energyNearest\", \"energyNextNearest\", \
\"muOrange\", \"muBlue\", \"inverseTemperature\".";

latticeGasBulkExchangeSweepCore::usage =
  "latticeGasBulkExchangeSweepCore[configuration, energyNearest, \
energyNextNearest, muOrange, muBlue, inverseTemperature, nAttempts] is the \
low-overhead fixed-argument form of latticeGasBulkExchangeSweep; returns \
{newConfiguration, statistics}.";

latticeGasBondSums::usage =
  "latticeGasBondSums[configuration] gives the exact integer bond sums \
{bondSumNearest, bondSumNextNearest} = {Sum s_i s_j over nearest-neighbor \
bonds, over diagonal bonds} (free boundaries), so the lattice energy is \
energyNearest*bondSumNearest + energyNextNearest*bondSumNextNearest.";

latticeGasTotalEnergy::usage =
  "latticeGasTotalEnergy[configuration, energyNearest, energyNextNearest] \
gives the total lattice energy.";

latticeGasSpeciesCounts::usage =
  "latticeGasSpeciesCounts[configuration] gives <|\"orange\" -> nOrange, \
\"blue\" -> nBlue, \"vacancy\" -> nVacancy|> for the real sites of the \
padded configuration.";

latticeGasConfigurationPlot::usage =
  "latticeGasConfigurationPlot[configuration, opts] shows the real sites of \
the padded configuration as an ArrayPlot (orange = +1, blue = -1, white = \
vacancy).  opts are passed to ArrayPlot.";

latticeGasOrderParameterMaps::usage =
  "latticeGasOrderParameterMaps[configuration, smoothingRadius] gives an \
Association of local order-parameter fields (real n x m matrices, values in \
[-1,1]), each a Gaussian-smoothed staggered average of the spin field:\n\
  \"composition\"      <s>            orange-rich (+) vs blue-rich (-);\n\
  \"occupancy\"        <|s|>          atoms (1) vs vacancies (0);\n\
  \"checkerboard\"     <(-1)^(i+j) s> checkerboard order, sign = registry;\n\
  \"stripeVertical\"   <(-1)^j s>     like-atom columns, sign = registry;\n\
  \"stripeHorizontal\" <(-1)^i s>     like-atom rows,    sign = registry.\n\
In a perfect single domain the matching map is +-1 and the others vanish; \
antiphase boundaries appear as sign changes (zero crossings).  \
smoothingRadius (default 8) should exceed the unit cell but be smaller than \
the domains; values within smoothingRadius of the boundary are damped by \
the zero padding.";

latticeGasOrderParameterPlot::usage =
  "latticeGasOrderParameterPlot[configuration, smoothingRadius, opts] shows \
the configuration beside its five order-parameter maps on a blue-white-\
orange diverging scale (-1 to +1).  opts are passed to ArrayPlot.";

latticeGasPhaseMap::usage =
  "latticeGasPhaseMap[configuration, smoothingRadius, orderThreshold] \
classifies every real site by its dominant local order parameter, returning \
an n x m matrix of integer codes: 0 vacancy region (occupancy < 1/2); \
1 orange-rich, 2 blue-rich; 3/4 checkerboard registries; 5/6 vertical-\
stripe registries; 7/8 horizontal-stripe registries; 9 disordered (no \
order parameter reaches orderThreshold, default 0.4).";

latticeGasPhaseMapPlot::usage =
  "latticeGasPhaseMapPlot[configuration, smoothingRadius, orderThreshold, \
opts] shows the phase map as a legended ArrayPlot: hue identifies the \
phase, light vs dark shade the registry (antiphase) variant, gray is \
disordered, white is vacancy.  opts are passed to ArrayPlot.";

latticeGasPhaseMapLegend::usage =
  "latticeGasPhaseMapLegend[] gives the phase-map SwatchLegend, for \
composite figures that share one legend among several \
latticeGasPhaseMapPlot panels (pass those plots \"Legend\" -> False).";

drawLattice::usage =
  "drawLattice[configuration] returns a graphics object of a square lattice with orange and blue disks corresponding
to +1 and -1 in the matrix configuration. Options: ImageSize->(Default 220)"

populationPie::usage =
  "populationPie[configuration] returns a graphics object of a pie chart indicating the populations of +1,-1, and 0: Options: ImageSize->(Default 100)"

Begin["`Private`"];

nearestNeighborKernel     = {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}};
nextNearestNeighborKernel = {{1, 0, 1}, {0, 0, 0}, {1, 0, 1}};

(* Requested counts are clipped rather than rejected, so interactive
   callers (sliders) can pass any combination: orange is capped first, blue
   fits into the remaining space, and at least one vacancy always survives
   (vacancy-mediated dynamics needs one to be ergodic). *)
latticeGasRandomConfiguration[nRows_Integer, nColumns_Integer,
   nOrangeAtoms_Integer, nBlueAtoms_Integer] :=
  Module[{nSites = nRows nColumns, nOrangeActual, nBlueActual},
   nOrangeActual = Clip[nOrangeAtoms, {0, nSites - 1}];
   nBlueActual = Clip[nBlueAtoms, {0, nSites - 1 - nOrangeActual}];
   ArrayPad[
    Partition[
     RandomSample[Join[
       ConstantArray[1, nOrangeActual],
       ConstantArray[-1, nBlueActual],
       ConstantArray[0, nSites - nOrangeActual - nBlueActual]]],
     nColumns], 1]];

(* ListConvolve with no overhang on the padded array returns exactly the
   neighbor sum at each real site; the zero frame supplies the free boundary. *)
latticeGasBondSums[configuration_] :=
  Module[{realSites = configuration[[2 ;; -2, 2 ;; -2]]},
   {Total[realSites ListConvolve[nearestNeighborKernel, configuration], 2]/2,
    Total[realSites ListConvolve[nextNearestNeighborKernel, configuration],
      2]/2}];

latticeGasTotalEnergy[configuration_, energyNearest_, energyNextNearest_] :=
  {energyNearest, energyNextNearest} . latticeGasBondSums[configuration];

latticeGasSpeciesCounts[configuration_] :=
  Module[{nRows, nColumns, nOrange, nBlue},
   {nRows, nColumns} = Dimensions[configuration] - 2;
   nOrange = Count[configuration, 1, {2}];
   nBlue = Count[configuration, -1, {2}];
   <|"orange" -> nOrange, "blue" -> nBlue,
    "vacancy" -> nRows nColumns - nOrange - nBlue|>];

latticeGasConfigurationPlot[configuration_, opts___] :=
  ArrayPlot[configuration[[2 ;; -2, 2 ;; -2]], opts,
   ColorRules -> {1 -> Orange, -1 -> RGBColor[0.25, 0.45, 0.85],
     0 -> White},
   Frame -> True, FrameTicks -> None, Mesh -> False];

(* Gaussian-smoothed staggered averages.  Multiplying by the sign pattern
   first maps each candidate phase onto a locally uniform field, which the
   smoothing then reads out; smoothing is done with zero padding, consistent
   with the free boundary. *)
smoothedField[field_, smoothingRadius_] :=
  ListConvolve[GaussianMatrix[smoothingRadius],
   ArrayPad[field, smoothingRadius]];

latticeGasOrderParameterMaps[configuration_, smoothingRadius_Integer: 8] :=
  Module[{realSites, nRows, nColumns, checkerboardSigns, columnSigns,
    rowSigns},
   realSites = N[configuration[[2 ;; -2, 2 ;; -2]]];
   {nRows, nColumns} = Dimensions[realSites];
   checkerboardSigns = Table[(-1.)^(i + j), {i, nRows}, {j, nColumns}];
   columnSigns = Table[(-1.)^j, {i, nRows}, {j, nColumns}];
   rowSigns = Table[(-1.)^i, {i, nRows}, {j, nColumns}];
   <|"composition" -> smoothedField[realSites, smoothingRadius],
    "occupancy" -> smoothedField[Abs[realSites], smoothingRadius],
    "checkerboard" ->
     smoothedField[checkerboardSigns realSites, smoothingRadius],
    "stripeVertical" ->
     smoothedField[columnSigns realSites, smoothingRadius],
    "stripeHorizontal" ->
     smoothedField[rowSigns realSites, smoothingRadius]|>];

divergingColorFunction =
  Blend[{RGBColor[0.25, 0.45, 0.85], White, Orange}, (# + 1)/2] &;

(* occupancy lives in [0,1] and is not a signed order parameter: use a
   grayscale (white = vacant, dark = occupied) so it cannot be confused
   with the orange/blue composition scale *)
occupancyColorFunction = GrayLevel[Abs[(1 - #)]^(1/2)] &;

style = Style[#,{FontSize->16, FontFamily->"Arial"}]&;

labelName = <|
"composition"-> style["Composition"],
"configuration"-> style["Microstructure"],
"occupancy"-> style["Vacancy Content"],
"checkerboard"-> style["Checkerboard Phase"],
"stripeVertical"->style["Stripes: Vertical"],
"stripeHorizontal"->style["Stripes: Horizontal"]
|>;


Options[drawLattice]={ImageSize-> 220};
drawLattice[configuration_, OptionsPattern[]]:= 
Module[
{orange =Position[configuration,1],blue =Position[configuration,-1] ,
rows,columns,
border = Transpose[{{0,0},Reverse@Dimensions[configuration]}]+ {{-1,1},{-1,1}}},
{rows, columns}= Dimensions[configuration];
orange ={ {0,1},{-1,0}} . # +{ 0,rows -1}&/@orange;
blue ={ {0,1},{-1,0}} . # +{ 0,rows -1}&/@blue;
Graphics[{{Orange,Disk[orange,0.45]},{Blue,Disk[blue,1/2]}, {FaceForm[],EdgeForm[LightGray],Rectangle[{1,-1},Reverse@Dimensions[configuration]-{0,2}]}},
ImageSize->OptionValue[ImageSize],(*Frame->True,*)
 PlotRange->border]
];

(* Counts only real sites (the padded frame would otherwise inflate the
   vacancy sector), and uses a fixed sector order so colors stay correct
   even when a species count is zero. *)
Options[populationPie]={ImageSize-> 100};
populationPie[configuration_, OptionsPattern[]] :=
 Module[{counts = latticeGasSpeciesCounts[configuration]},
  PieChart[
   {Labeled[counts["blue"], "blue"], Labeled[counts["vacancy"], "vacancy"],
    Labeled[counts["orange"], "orange"]},
   ChartStyle -> {RGBColor[0.25, 0.45, 0.85], White, Orange},
   ImageSize -> OptionValue[ImageSize], PerformanceGoal -> "Speed"]]


latticeGasOrderParameterPlot[configuration_,
   smoothingRadius_Integer: 8, opts___] :=
  Module[{maps = latticeGasOrderParameterMaps[configuration,
      smoothingRadius]},
   Grid[Partition[
     Join[
      {Labeled[
      (*latticeGasConfigurationPlot[configuration, opts,  ImageSize -> 220]*)
      drawLattice[configuration], 
         labelName["configuration"], Top]},
      Table[
       Labeled[ArrayPlot[maps[mapName], opts,
         ColorFunction -> If[mapName === "occupancy",
           occupancyColorFunction, divergingColorFunction],
         ColorFunctionScaling -> False,
         PlotRange -> If[mapName === "occupancy", {0, 1}, {-1, 1}],
         Frame -> True, FrameTicks -> None, ImageSize -> 220],
        labelName[mapName],Top],
       {mapName, {"composition", "occupancy", "checkerboard",
         "stripeVertical", "stripeHorizontal"}}]],
     3], Spacings -> {1, 1}]];

(* phase map: hue = phase, light/dark shade = registry variant *)
phaseMapColorRules = {
   0 -> White,                            (* vacancy region *)
   1 -> Orange,                           (* orange-rich *)
   2 -> RGBColor[0.25, 0.45, 0.85],       (* blue-rich *)
   3 -> RGBColor[0.13, 0.55, 0.13],       (* checkerboard, registry A *)
   4 -> RGBColor[0.62, 0.85, 0.50],       (* checkerboard, registry B *)
   5 -> RGBColor[0.72, 0.15, 0.15],       (* vertical stripes, registry A *)
   6 -> RGBColor[0.98, 0.62, 0.60],       (* vertical stripes, registry B *)
   7 -> RGBColor[0.42, 0.22, 0.65],       (* horizontal stripes, registry A *)
   8 -> RGBColor[0.80, 0.68, 0.90],       (* horizontal stripes, registry B *)
   9 -> GrayLevel[0.88]};                 (* disordered *)

phaseMapLabels = {"orange-rich", "blue-rich",
   "checkerboard (registry A)", "checkerboard (registry B)",
   "vertical stripes (registry A)", "vertical stripes (registry B)",
   "horizontal stripes (registry A)", "horizontal stripes (registry B)",
   "disordered", "vacancy"};

latticeGasPhaseMapLegend[] :=
  SwatchLegend[Lookup[Association[phaseMapColorRules],
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 0}], phaseMapLabels];

latticeGasPhaseMap[configuration_, smoothingRadius_Integer: 8,
   Optional[orderThreshold_?NumericQ, 0.4]] :=
  Module[{maps, composition, occupancy, checkerboard, stripeVertical,
    stripeHorizontal, nRows, nColumns},
   maps = latticeGasOrderParameterMaps[configuration, smoothingRadius];
   composition = maps["composition"]; occupancy = maps["occupancy"];
   checkerboard = maps["checkerboard"];
   stripeVertical = maps["stripeVertical"];
   stripeHorizontal = maps["stripeHorizontal"];
   {nRows, nColumns} = Dimensions[composition];
   Table[
    Module[{values, dominantIndex},
     If[occupancy[[i, j]] < 0.5, 0,
      values = {composition[[i, j]], checkerboard[[i, j]],
        stripeVertical[[i, j]], stripeHorizontal[[i, j]]};
      dominantIndex = First[Ordering[Abs[values], -1]];
      If[Abs[values[[dominantIndex]]] < orderThreshold, 9,
       2 dominantIndex - If[values[[dominantIndex]] > 0, 1, 0]]]],
    {i, nRows}, {j, nColumns}]];

latticeGasPhaseMapPlot[configuration_, smoothingRadius_Integer: 8,
   Optional[orderThreshold_?NumericQ, 0.4], opts___] :=
  Module[{includeLegend, arrayPlotOptions, plot},
   includeLegend = ! MemberQ[{opts}, "Legend" -> False];
   arrayPlotOptions = DeleteCases[{opts}, "Legend" -> _];
   plot = ArrayPlot[
     latticeGasPhaseMap[configuration, smoothingRadius, orderThreshold],
     Sequence @@ arrayPlotOptions,
     ColorRules -> phaseMapColorRules,
     Frame -> True, FrameTicks -> None];
   If[includeLegend, Legended[plot, latticeGasPhaseMapLegend[]], plot]];

$latticeGasKernel = With[
   {rowOffsetTable = {-1, 1, 0, 0, -1, -1, 1, 1},
    columnOffsetTable = {0, 0, -1, 1, -1, 1, -1, 1}},
   Compile[{{initialConfiguration, _Integer, 2},
     {energyNearest, _Real}, {energyNextNearest, _Real},
     {muOrange, _Real}, {muBlue, _Real},
     {inverseTemperature, _Real}, {exchangeProbability, _Real},
     {vacancySwitchFrequency, _Real}, {occupiedSwitchFrequency, _Real},
     {nAttempts, _Integer}},
    Module[{configuration = initialConfiguration,
      nRows, nColumns, nSites, surfaceSiteCount,
      vacancyRowList, vacancyColumnList,
      occupiedRowList, occupiedColumnList, listSlot,
      vacancyCount = 0, occupiedCount = 0,
      i = 1, j = 1, neighborRow, neighborColumn,
      surfaceIndex, offsetIndex, rowOffset, columnOffset,
      siteValue, neighborValue, proposedSpecies, currentSiteValue,
      siteEnthalpy, neighborEnthalpy,
      bondCoupling, deltaEnergy, muProposedSpecies,
      vacancyWeight, occupiedWeight, totalWeight,
      chosenSlot, atomSlot, vacancySlot, lastRow, lastColumn,
      attemptedVacancySwitches = 0, acceptedVacancySwitches = 0,
      attemptedOccupiedSwitches = 0, acceptedOccupiedSwitches = 0,
      attemptedExchanges = 0, acceptedInsertions = 0, acceptedRemovals = 0,
      statisticsRow1, statisticsRow2},
     nRows = Length[configuration] - 2;
     nColumns = Length[First[configuration]] - 2;
     nSites = nRows nColumns;
     surfaceSiteCount = 2 nRows + 2 nColumns - 4;
     (* ---- build the vacancy and atom position lists (one O(N M) scan);
        listSlot[[i,j]] records where site (i,j) sits in its class list ---- *)
     vacancyRowList = Table[0, {nSites}];
     vacancyColumnList = Table[0, {nSites}];
     occupiedRowList = Table[0, {nSites}];
     occupiedColumnList = Table[0, {nSites}];
     listSlot = Table[0, {nRows + 2}, {nColumns + 2}];
     Do[
      If[configuration[[i, j]] == 0,
       vacancyCount++;
       vacancyRowList[[vacancyCount]] = i;
       vacancyColumnList[[vacancyCount]] = j;
       listSlot[[i, j]] = vacancyCount,
       occupiedCount++;
       occupiedRowList[[occupiedCount]] = i;
       occupiedColumnList[[occupiedCount]] = j;
       listSlot[[i, j]] = occupiedCount],
      {i, 2, nRows + 1}, {j, 2, nColumns + 1}];
     Do[
      If[RandomReal[] < exchangeProbability,
       (* ---- surface exchange with the gas reservoir ---- *)
       attemptedExchanges++;
       surfaceIndex = RandomInteger[{1, surfaceSiteCount}];
       Which[
        surfaceIndex <= nColumns,
        i = 2; j = surfaceIndex + 1,
        surfaceIndex <= 2 nColumns,
        i = nRows + 1; j = surfaceIndex - nColumns + 1,
        surfaceIndex <= 2 nColumns + nRows - 2,
        i = surfaceIndex - 2 nColumns + 2; j = 2,
        True,
        i = surfaceIndex - 2 nColumns - nRows + 4; j = nColumns + 1];
       proposedSpecies = 2 RandomInteger[{0, 1}] - 1;
       muProposedSpecies = If[proposedSpecies == 1, muOrange, muBlue];
       currentSiteValue = configuration[[i, j]];
       If[currentSiteValue == proposedSpecies || currentSiteValue == 0,
        siteEnthalpy =
         energyNearest (configuration[[i - 1, j]] +
            configuration[[i + 1, j]] + configuration[[i, j - 1]] +
            configuration[[i, j + 1]]) +
         energyNextNearest (configuration[[i - 1, j - 1]] +
            configuration[[i - 1, j + 1]] + configuration[[i + 1, j - 1]] +
            configuration[[i + 1, j + 1]]);
        If[currentSiteValue == proposedSpecies,
         (* removal to the gas *)
         deltaEnergy = muProposedSpecies - proposedSpecies siteEnthalpy;
         If[deltaEnergy <= 0. ||
           RandomReal[] < Exp[-inverseTemperature deltaEnergy],
          configuration[[i, j]] = 0; acceptedRemovals++;
          (* atom list: overwrite this atom's slot with the last entry *)
          atomSlot = listSlot[[i, j]];
          lastRow = occupiedRowList[[occupiedCount]];
          lastColumn = occupiedColumnList[[occupiedCount]];
          occupiedRowList[[atomSlot]] = lastRow;
          occupiedColumnList[[atomSlot]] = lastColumn;
          listSlot[[lastRow, lastColumn]] = atomSlot;
          occupiedCount--;
          (* vacancy list: append the new vacancy *)
          vacancyCount++;
          vacancyRowList[[vacancyCount]] = i;
          vacancyColumnList[[vacancyCount]] = j;
          listSlot[[i, j]] = vacancyCount],
         (* insertion from the gas *)
         deltaEnergy = proposedSpecies siteEnthalpy - muProposedSpecies;
         If[deltaEnergy <= 0. ||
           RandomReal[] < Exp[-inverseTemperature deltaEnergy],
          configuration[[i, j]] = proposedSpecies; acceptedInsertions++;
          (* vacancy list: overwrite this vacancy's slot with the last entry *)
          vacancySlot = listSlot[[i, j]];
          lastRow = vacancyRowList[[vacancyCount]];
          lastColumn = vacancyColumnList[[vacancyCount]];
          vacancyRowList[[vacancySlot]] = lastRow;
          vacancyColumnList[[vacancySlot]] = lastColumn;
          listSlot[[lastRow, lastColumn]] = vacancySlot;
          vacancyCount--;
          (* atom list: append the new atom *)
          occupiedCount++;
          occupiedRowList[[occupiedCount]] = i;
          occupiedColumnList[[occupiedCount]] = j;
          listSlot[[i, j]] = occupiedCount]]],
       (* ---- switch attempt: choose class with probability proportional to
          frequency * population, so per-site rates are composition-independent *)
       vacancyWeight = vacancySwitchFrequency vacancyCount;
       occupiedWeight = occupiedSwitchFrequency occupiedCount;
       totalWeight = vacancyWeight + occupiedWeight;
       If[totalWeight > 0.,
        If[RandomReal[] totalWeight < vacancyWeight,
         (* ---- vacancy switch: vacancy <-> neighboring atom ---- *)
         attemptedVacancySwitches++;
         chosenSlot = RandomInteger[{1, vacancyCount}];
         i = vacancyRowList[[chosenSlot]];
         j = vacancyColumnList[[chosenSlot]];
         offsetIndex = RandomInteger[{1, 8}];
         rowOffset = rowOffsetTable[[offsetIndex]];
         columnOffset = columnOffsetTable[[offsetIndex]];
         neighborRow = i + rowOffset; neighborColumn = j + columnOffset;
         If[neighborRow >= 2 && neighborRow <= nRows + 1 &&
           neighborColumn >= 2 && neighborColumn <= nColumns + 1,
          neighborValue = configuration[[neighborRow, neighborColumn]];
          If[neighborValue != 0,
           siteValue = 0;
           siteEnthalpy =
            energyNearest (configuration[[i - 1, j]] +
               configuration[[i + 1, j]] + configuration[[i, j - 1]] +
               configuration[[i, j + 1]]) +
            energyNextNearest (configuration[[i - 1, j - 1]] +
               configuration[[i - 1, j + 1]] +
               configuration[[i + 1, j - 1]] +
               configuration[[i + 1, j + 1]]);
           neighborEnthalpy =
            energyNearest (configuration[[neighborRow - 1,
                neighborColumn]] +
               configuration[[neighborRow + 1, neighborColumn]] +
               configuration[[neighborRow, neighborColumn - 1]] +
               configuration[[neighborRow, neighborColumn + 1]]) +
            energyNextNearest (configuration[[neighborRow - 1,
                neighborColumn - 1]] +
               configuration[[neighborRow - 1, neighborColumn + 1]] +
               configuration[[neighborRow + 1, neighborColumn - 1]] +
               configuration[[neighborRow + 1, neighborColumn + 1]]);
           bondCoupling = If[rowOffset columnOffset == 0,
             energyNearest, energyNextNearest];
           deltaEnergy = (siteValue - neighborValue) *
             ((neighborEnthalpy - bondCoupling siteValue) -
              (siteEnthalpy - bondCoupling neighborValue));
           If[deltaEnergy <= 0. ||
             RandomReal[] < Exp[-inverseTemperature deltaEnergy],
            configuration[[i, j]] = neighborValue;
            configuration[[neighborRow, neighborColumn]] = 0;
            acceptedVacancySwitches++;
            (* the vacancy and the atom trade places in the lists *)
            atomSlot = listSlot[[neighborRow, neighborColumn]];
            occupiedRowList[[atomSlot]] = i;
            occupiedColumnList[[atomSlot]] = j;
            vacancyRowList[[chosenSlot]] = neighborRow;
            vacancyColumnList[[chosenSlot]] = neighborColumn;
            listSlot[[i, j]] = atomSlot;
            listSlot[[neighborRow, neighborColumn]] = chosenSlot]]],
         (* ---- occupied switch: atom <-> neighboring atom of the other
            species (a vacancy neighbor is NOT eligible: that channel
            belongs exclusively to the vacancy class) ---- *)
         attemptedOccupiedSwitches++;
         chosenSlot = RandomInteger[{1, occupiedCount}];
         i = occupiedRowList[[chosenSlot]];
         j = occupiedColumnList[[chosenSlot]];
         offsetIndex = RandomInteger[{1, 8}];
         rowOffset = rowOffsetTable[[offsetIndex]];
         columnOffset = columnOffsetTable[[offsetIndex]];
         neighborRow = i + rowOffset; neighborColumn = j + columnOffset;
         If[neighborRow >= 2 && neighborRow <= nRows + 1 &&
           neighborColumn >= 2 && neighborColumn <= nColumns + 1,
          siteValue = configuration[[i, j]];
          neighborValue = configuration[[neighborRow, neighborColumn]];
          If[neighborValue != 0 && neighborValue != siteValue,
           siteEnthalpy =
            energyNearest (configuration[[i - 1, j]] +
               configuration[[i + 1, j]] + configuration[[i, j - 1]] +
               configuration[[i, j + 1]]) +
            energyNextNearest (configuration[[i - 1, j - 1]] +
               configuration[[i - 1, j + 1]] +
               configuration[[i + 1, j - 1]] +
               configuration[[i + 1, j + 1]]);
           neighborEnthalpy =
            energyNearest (configuration[[neighborRow - 1,
                neighborColumn]] +
               configuration[[neighborRow + 1, neighborColumn]] +
               configuration[[neighborRow, neighborColumn - 1]] +
               configuration[[neighborRow, neighborColumn + 1]]) +
            energyNextNearest (configuration[[neighborRow - 1,
                neighborColumn - 1]] +
               configuration[[neighborRow - 1, neighborColumn + 1]] +
               configuration[[neighborRow + 1, neighborColumn - 1]] +
               configuration[[neighborRow + 1, neighborColumn + 1]]);
           bondCoupling = If[rowOffset columnOffset == 0,
             energyNearest, energyNextNearest];
           deltaEnergy = (siteValue - neighborValue) *
             ((neighborEnthalpy - bondCoupling siteValue) -
              (siteEnthalpy - bondCoupling neighborValue));
           If[deltaEnergy <= 0. ||
             RandomReal[] < Exp[-inverseTemperature deltaEnergy],
            configuration[[i, j]] = neighborValue;
            configuration[[neighborRow, neighborColumn]] = siteValue;
            acceptedOccupiedSwitches++
            (* both sites remain occupied: the lists are unchanged *)]]]]]],
      {nAttempts}];
     statisticsRow1 = Table[0, {nColumns + 2}];
     statisticsRow1[[1]] = attemptedVacancySwitches;
     statisticsRow1[[2]] = acceptedVacancySwitches;
     statisticsRow1[[3]] = attemptedOccupiedSwitches;
     statisticsRow1[[4]] = acceptedOccupiedSwitches;
     statisticsRow1[[5]] = attemptedExchanges;
     statisticsRow2 = Table[0, {nColumns + 2}];
     statisticsRow2[[1]] = acceptedInsertions;
     statisticsRow2[[2]] = acceptedRemovals;
     Join[configuration, {statisticsRow1}, {statisticsRow2}]],
    CompilationTarget -> "C", RuntimeOptions -> "Speed"]];

statisticsKeys = {"attemptedVacancySwitches", "acceptedVacancySwitches",
   "attemptedOccupiedSwitches", "acceptedOccupiedSwitches",
   "attemptedExchanges", "acceptedInsertions", "acceptedRemovals"};

latticeGasMetropolisSweepCore[configuration_, energyNearest_,
   energyNextNearest_, muOrange_, muBlue_, inverseTemperature_,
   exchangeProbability_, vacancySwitchFrequency_,
   occupiedSwitchFrequency_, nAttempts_] :=
  Module[{result = $latticeGasKernel[configuration, energyNearest,
      energyNextNearest, muOrange, muBlue, inverseTemperature,
      exchangeProbability, vacancySwitchFrequency,
      occupiedSwitchFrequency, nAttempts]},
   {result[[;; -3]],
    AssociationThread[statisticsKeys,
     Join[result[[-2, 1 ;; 5]], result[[-1, 1 ;; 2]]]]}];

Options[latticeGasMetropolisSweep] = {"energyNearest" -> 1.,
   "energyNextNearest" -> 0., "muOrange" -> 0., "muBlue" -> 0.,
   "inverseTemperature" -> 1., "exchangeProbability" -> 0.,
   "vacancySwitchFrequency" -> 1., "occupiedSwitchFrequency" -> 1.};

latticeGasMetropolisSweep[configuration_, nAttempts_Integer,
   OptionsPattern[]] :=
  latticeGasMetropolisSweepCore[configuration,
   OptionValue["energyNearest"],
   OptionValue["energyNextNearest"],
   OptionValue["muOrange"],
   OptionValue["muBlue"],
   OptionValue["inverseTemperature"],
   OptionValue["exchangeProbability"],
   OptionValue["vacancySwitchFrequency"],
   OptionValue["occupiedSwitchFrequency"], nAttempts];

(* ---- bulk-exchange (phase-diagram) tool ---- *)

$latticeGasBulkKernel =
  Compile[{{initialConfiguration, _Integer, 2},
    {energyNearest, _Real}, {energyNextNearest, _Real},
    {muOrange, _Real}, {muBlue, _Real},
    {inverseTemperature, _Real}, {nAttempts, _Integer}},
   Module[{configuration = initialConfiguration,
     nRows, nColumns, i, j, currentSiteValue, proposedValue,
     siteEnthalpy, deltaEnergy, muCurrent, muProposed,
     attemptedBulkExchanges = 0, acceptedToOrange = 0,
     acceptedToBlue = 0, acceptedToVacancy = 0, statisticsRow},
    nRows = Length[configuration] - 2;
    nColumns = Length[First[configuration]] - 2;
    Do[
     attemptedBulkExchanges++;
     i = RandomInteger[{2, nRows + 1}];
     j = RandomInteger[{2, nColumns + 1}];
     currentSiteValue = configuration[[i, j]];
     (* propose one of the OTHER two states, each with probability 1/2:
        a symmetric proposal, so plain Metropolis acceptance is exact *)
     proposedValue =
      Mod[currentSiteValue + 2 + RandomInteger[{0, 1}], 3] - 1;
     muCurrent = If[currentSiteValue == 1, muOrange,
       If[currentSiteValue == -1, muBlue, 0.]];
     muProposed = If[proposedValue == 1, muOrange,
       If[proposedValue == -1, muBlue, 0.]];
     siteEnthalpy =
      energyNearest (configuration[[i - 1, j]] +
         configuration[[i + 1, j]] + configuration[[i, j - 1]] +
         configuration[[i, j + 1]]) +
      energyNextNearest (configuration[[i - 1, j - 1]] +
         configuration[[i - 1, j + 1]] + configuration[[i + 1, j - 1]] +
         configuration[[i + 1, j + 1]]);
     deltaEnergy = (proposedValue - currentSiteValue) siteEnthalpy -
       (muProposed - muCurrent);
     If[deltaEnergy <= 0. ||
       RandomReal[] < Exp[-inverseTemperature deltaEnergy],
      configuration[[i, j]] = proposedValue;
      Which[
       proposedValue == 1, acceptedToOrange++,
       proposedValue == -1, acceptedToBlue++,
       True, acceptedToVacancy++]],
     {nAttempts}];
    statisticsRow = Table[0, {nColumns + 2}];
    statisticsRow[[1]] = attemptedBulkExchanges;
    statisticsRow[[2]] = acceptedToOrange + acceptedToBlue +
      acceptedToVacancy;
    statisticsRow[[3]] = acceptedToOrange;
    statisticsRow[[4]] = acceptedToBlue;
    statisticsRow[[5]] = acceptedToVacancy;
    Join[configuration, {statisticsRow}]],
   CompilationTarget -> "C", RuntimeOptions -> "Speed"];

bulkStatisticsKeys = {"attemptedBulkExchanges", "acceptedBulkExchanges",
   "acceptedToOrange", "acceptedToBlue", "acceptedToVacancy"};

latticeGasBulkExchangeSweepCore[configuration_, energyNearest_,
   energyNextNearest_, muOrange_, muBlue_, inverseTemperature_,
   nAttempts_] :=
  Module[{result = $latticeGasBulkKernel[configuration, energyNearest,
      energyNextNearest, muOrange, muBlue, inverseTemperature,
      nAttempts]},
   {result[[;; -2]],
    AssociationThread[bulkStatisticsKeys, result[[-1, 1 ;; 5]]]}];

Options[latticeGasBulkExchangeSweep] = {"energyNearest" -> 1.,
   "energyNextNearest" -> 0., "muOrange" -> 0., "muBlue" -> 0.,
   "inverseTemperature" -> 1., "exchangeProbability"->None};

latticeGasBulkExchangeSweep[configuration_, nAttempts_Integer,
   OptionsPattern[]] :=
  latticeGasBulkExchangeSweepCore[configuration,
   OptionValue["energyNearest"],
   OptionValue["energyNextNearest"],
   OptionValue["muOrange"],
   OptionValue["muBlue"],
   OptionValue["inverseTemperature"], nAttempts];

End[];
EndPackage[];
