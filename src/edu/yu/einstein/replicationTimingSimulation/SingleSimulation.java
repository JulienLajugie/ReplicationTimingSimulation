/*******************************************************************************
 *     GenPlay, Einstein Genome Analyzer
 *     Copyright (C) 2009, 2011 Albert Einstein College of Medicine
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *     Authors:	Julien Lajugie <julien.lajugie@einstein.yu.edu>
 *     			Nicolas Fourel <nicolas.fourel@einstein.yu.edu>
 *     Website: <http://genplay.einstein.yu.edu>
 *******************************************************************************/
package edu.yu.einstein.replicationTimingSimulation;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import edu.yu.einstein.genplay.core.operation.Operation;
import edu.yu.einstein.genplay.core.operation.SCWList.MCWLOInvertMask;
import edu.yu.einstein.genplay.core.operation.SCWList.SCWLOConvertIntoBinList;
import edu.yu.einstein.genplay.core.operation.SCWList.SCWLOConvertIntoGeneList;
import edu.yu.einstein.genplay.core.operation.SCWList.SCWLOConvertIntoSimpleSCWList;
import edu.yu.einstein.genplay.core.operation.SCWList.SCWLOFilterThreshold;
import edu.yu.einstein.genplay.core.operation.SCWList.SCWLOOperationWithConstant;
import edu.yu.einstein.genplay.core.operation.SCWList.SCWLOTwoLayers;
import edu.yu.einstein.genplay.core.operation.binList.BLOFindIslands;
import edu.yu.einstein.genplay.core.operation.binList.BLOGauss;
import edu.yu.einstein.genplay.core.operation.binList.BLOTwoLayers;
import edu.yu.einstein.genplay.core.operation.geneList.GLOMergeGeneLists;
import edu.yu.einstein.genplay.core.operation.geneList.GLOScoreFromSCWList;
import edu.yu.einstein.genplay.dataStructure.enums.GeneScoreType;
import edu.yu.einstein.genplay.dataStructure.enums.IslandResultType;
import edu.yu.einstein.genplay.dataStructure.enums.OperationWithConstant;
import edu.yu.einstein.genplay.dataStructure.enums.SCWListType;
import edu.yu.einstein.genplay.dataStructure.enums.ScoreOperation;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.SCWList.SCWList;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.SCWList.binList.BinList;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.geneList.GeneList;
import edu.yu.einstein.genplay.dataStructure.list.listView.ListView;
import edu.yu.einstein.genplay.dataStructure.scoredChromosomeWindow.ScoredChromosomeWindow;


/**
 * Runs a simulation for a given island size and a given
 * @author Julien Lajugie
 */
public class SingleSimulation implements Operation<SimulationResult> {

	private final static int		ISLAND_DISTANCE = 4000000; 	// space between 2 island starts position

	private final static float		Q_VALUE_CUTOFF = 0.05f;		// cutoff for the qValue
	private final static boolean	USE_ISLAND_FINDER = true;	// true to use the island finder to define the island
	/**
	 * Prints the specified {@link ListView} in a temporary file with the specified name prefix
	 * @param data
	 * @param prefix
	 * @throws IOException
	 */
	public static void printSCWInTmpFile(ListView<? extends ScoredChromosomeWindow> data, String prefix) throws IOException {
		File file = File.createTempFile(prefix, ".bed");
		System.out.println("tmp file: " + file.getPath());
		//Write the data
		PrintWriter dout = new PrintWriter(file);
		for (int i = 0; i < data.size(); i++) {
			ScoredChromosomeWindow scw = data.get(i);
			dout.println("chr1\t" + (scw.getStart() - 1)+ "\t" + (scw.getStop() - 1) + "\t-\t" + scw.getScore());
		}
		dout.close();
	}

	private final int 		islandSize;					// size of the islands to use in the simulation
	private final double 	percentageReadToAdd;		// percentage of reads to add in the S phase in the islands
	private final SCWList 	sList;						// s phase data
	private final SCWList 	g1List;						// g1 phase data


	/**
	 * Creates an instance of {@link SingleSimulation}
	 * @param islandSize size of the islands to use in the simulation
	 * @param percentageReadToAdd percentage of reads to add in the S phase in the islands
	 * @param sList s phase data
	 * @param g1List g1 phase data
	 */
	public SingleSimulation(int islandSize,
			double percentageReadToAdd,
			SCWList sList,
			SCWList g1List) {
		this.islandSize = islandSize;
		this.percentageReadToAdd = percentageReadToAdd;
		this.sList = sList;
		this.g1List = g1List;
	}

	@Override
	public SimulationResult compute() throws Exception {
		System.out.println("SingleSimulation.compute() - 1");
		// 1 - generate control lists
		SCWList[] resampledList = new ResampleLayers(sList, g1List, 0).compute();
		SCWList controlS = resampledList[0];
		SCWList controlG1 = resampledList[1];

		// 2 - generate sample lists

		// 2a - generate list with reads added
		System.out.println("SingleSimulation.compute() - 2a");
		resampledList = new ResampleLayers(sList, g1List, percentageReadToAdd).compute();
		SCWList resampledSReadAdded = resampledList[0];
		SCWList resampledG1ReadAdded = resampledList[1];

		// 2b - generate list with no reads added
		System.out.println("SingleSimulation.compute() - 2b");
		resampledList = new ResampleLayers(sList, g1List, 0).compute();
		SCWList resampledSNoReadAdded = resampledList[0];
		SCWList resampledG1NoReadAdded = resampledList[1];

		// 2c - generate islands mask list and its inverse
		System.out.println("SingleSimulation.compute() - 2c");
		SCWList islandMask = new GenerateIslands(ISLAND_DISTANCE, islandSize).compute();
		SCWList islandInvertedMask = new MCWLOInvertMask(islandMask).compute();
		//printSCWInTmpFile(islandMask.get(0), "inislands_" + islandSize + "_" + percentageReadToAdd + "_");

		// 2d - sum the island with reads added with the baseline
		System.out.println("SingleSimulation.compute() - 2d");
		resampledSReadAdded = new SCWLOTwoLayers(resampledSReadAdded, islandMask, ScoreOperation.MULTIPLICATION).compute();
		resampledSNoReadAdded = new SCWLOTwoLayers(resampledSNoReadAdded, islandInvertedMask, ScoreOperation.MULTIPLICATION).compute();
		SCWList resampledS = new SCWLOTwoLayers(resampledSReadAdded, resampledSNoReadAdded, ScoreOperation.ADDITION).compute();

		resampledG1ReadAdded = new SCWLOTwoLayers(resampledG1ReadAdded, islandMask, ScoreOperation.MULTIPLICATION).compute();
		resampledG1NoReadAdded = new SCWLOTwoLayers(resampledG1NoReadAdded, islandInvertedMask, ScoreOperation.MULTIPLICATION).compute();
		SCWList resampledG1 = new SCWLOTwoLayers(resampledG1ReadAdded, resampledG1NoReadAdded, ScoreOperation.ADDITION).compute();

		// 3 - convert into binlist
		System.out.println("SingleSimulation.compute() - 3");
		BinList binnedControlS = new SCWLOConvertIntoBinList(controlS, 500, ScoreOperation.ADDITION).compute();
		BinList binnedControlG1 = new SCWLOConvertIntoBinList(controlG1, 500, ScoreOperation.ADDITION).compute();
		BinList binnedResampledS = new SCWLOConvertIntoBinList(resampledS, 500, ScoreOperation.ADDITION).compute();
		BinList binnedResampledG1 = new SCWLOConvertIntoBinList(resampledG1, 500, ScoreOperation.ADDITION).compute();

		// 4 - gauss binlists
		System.out.println("SingleSimulation.compute() - 4");
		binnedControlS = new BLOGauss(binnedControlS, 400000, false).compute();
		binnedControlG1 = new BLOGauss(binnedControlG1, 400000, false).compute();
		binnedResampledS = new BLOGauss(binnedResampledS, 400000, false).compute();
		binnedResampledG1 = new BLOGauss(binnedResampledG1, 400000, false).compute();

		// 5 - compute S / G1 ratios
		System.out.println("SingleSimulation.compute() - 5");
		BinList controlSG1 = (BinList) new BLOTwoLayers(binnedControlS, binnedControlG1, ScoreOperation.DIVISION).compute();
		BinList sampleSG1 = (BinList) new BLOTwoLayers(binnedResampledS, binnedResampledG1, ScoreOperation.DIVISION).compute();

		// 6 - remove windows that are null in one of the 2 lists
		System.out.println("SingleSimulation.compute() - 6");
		BinList maskControl = (BinList) new SCWLOOperationWithConstant(controlSG1, OperationWithConstant.UNIQUE_SCORE, 1, false).compute();
		BinList maskSample = (BinList) new SCWLOOperationWithConstant(sampleSG1, OperationWithConstant.UNIQUE_SCORE, 1, false).compute();
		controlSG1 = (BinList) new BLOTwoLayers(maskSample, controlSG1, ScoreOperation.MULTIPLICATION).compute();
		sampleSG1 = (BinList) new BLOTwoLayers(maskControl, sampleSG1, ScoreOperation.MULTIPLICATION).compute();
		//printSCWInTmpFile(controlSG1.get(0), "controlSG1_" + islandSize + "_" + percentageReadToAdd + "_");
		//printSCWInTmpFile(sampleSG1.get(0), "sampleSG1_" + islandSize + "_" + percentageReadToAdd + "_");

		// 7 - compute sample - control difference
		System.out.println("SingleSimulation.compute() - 7");
		BinList sampleCtrlDifference = (BinList) new BLOTwoLayers(sampleSG1, controlSG1, ScoreOperation.SUBTRACTION).compute();
		//printSCWInTmpFile(sampleCtrlDifference.get(0), "difference_" + islandSize + "_" + percentageReadToAdd + "_");

		// 8 - call islands
		System.out.println("SingleSimulation.compute() - 8");
		GeneList islands;
		if (USE_ISLAND_FINDER) {
			BinList positiveSampleCtrlDifference = (BinList) new SCWLOFilterThreshold(sampleCtrlDifference, 0, Float.POSITIVE_INFINITY, false).compute();
			BinList negativesSampleCtrlDifference = (BinList) new SCWLOOperationWithConstant(sampleCtrlDifference, OperationWithConstant.MULTIPLICATION, -1f, false).compute();
			negativesSampleCtrlDifference = (BinList) new SCWLOFilterThreshold(negativesSampleCtrlDifference, 0, Float.POSITIVE_INFINITY, false).compute();
			GeneList positiveIslands = findIslandUsingIslandFinder(positiveSampleCtrlDifference);
			GeneList negativeIslands = findIslandUsingIslandFinder(negativesSampleCtrlDifference);
			islands = new GLOMergeGeneLists(positiveIslands, negativeIslands).compute();

		} else {
			islands = new FindIslands(sampleCtrlDifference).compute();
		}

		// 9 - score islands
		System.out.println("SingleSimulation.compute() - 9");
		GeneList controlIslandsS = new GLOScoreFromSCWList(islands, controlS, GeneScoreType.BASE_COVERAGE_SUM).compute();
		GeneList controlIslandsG1 = new GLOScoreFromSCWList(islands, controlG1, GeneScoreType.BASE_COVERAGE_SUM).compute();
		GeneList sampleIslandsS = new GLOScoreFromSCWList(islands, resampledS, GeneScoreType.BASE_COVERAGE_SUM).compute();
		GeneList sampleIslandsG1 = new GLOScoreFromSCWList(islands, resampledG1, GeneScoreType.BASE_COVERAGE_SUM).compute();

		// 10 - compute fisher exact test and retrieve qvalues
		System.out.println("SingleSimulation.compute() - 10");
		SCWList islandsQValues = new ComputeQValues(controlIslandsS, controlIslandsG1, sampleIslandsS, sampleIslandsG1).compute();
		//printSCWInTmpFile(islandsQValues.get(0), "islands_before_filter" + islandSize + "_" + percentageReadToAdd + "_");

		// 11 - filter islands with qvalue under 0.05
		System.out.println("SingleSimulation.compute() - 11");
		SCWList filteredIslands = new SCWLOFilterThreshold(islandsQValues, 0, Q_VALUE_CUTOFF, false).compute();
		filteredIslands = new SCWLOConvertIntoSimpleSCWList(filteredIslands, SCWListType.MASK).compute();
		//printSCWInTmpFile(filteredIslands.get(0), "islands_after_filter" + islandSize + "_" + percentageReadToAdd + "_");

		// 12 - compute false positive and false negatives
		System.out.println("SingleSimulation.compute() - 12");
		SimulationResult simulationResult = new ComputeSimulationResult(islandSize, percentageReadToAdd, islandMask, filteredIslands).compute();

		return simulationResult;
	}


	/**
	 * Find the islands in the input genelist using the genplay island finder algorithm
	 * @param input
	 * @return
	 * @throws Exception
	 */
	private GeneList findIslandUsingIslandFinder(BinList input) throws Exception {
		BLOFindIslands bloFindIslands = new BLOFindIslands(input);
		bloFindIslands.getIsland().setWindowMinValue(0.1);
		bloFindIslands.getIsland().setGap(100);
		bloFindIslands.getIsland().setIslandMinScore(0);
		bloFindIslands.getIsland().setIslandMinLength(100);
		IslandResultType[] resType = {IslandResultType.IFSCORE};
		bloFindIslands.setList(resType);
		BinList resBinList = bloFindIslands.compute()[0];
		SCWList resMaskList = new SCWLOConvertIntoSimpleSCWList(resBinList, SCWListType.MASK).compute();
		GeneList resGeneList = new SCWLOConvertIntoGeneList(resMaskList).compute();
		return resGeneList;
	}


	@Override
	public String getDescription() {
		return "Operation: Run Simulation";
	}


	@Override
	public String getProcessingDescription() {
		return "Running Simmulation";
	}


	@Override
	public int getStepCount() {
		return 0;
	}


	@Override
	public void stop() {
		// stop operation not implemented
	}
}
