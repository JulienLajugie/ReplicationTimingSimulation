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

	private final static int		ISLAND_DISTANCE 	= 4000000; 	// space between 2 island starts position
	private final static float		Q_VALUE_CUTOFF 		= 0.05f;	// cutoff for the qValue
	private final static boolean	USE_ISLAND_FINDER 	= true;		// true to use the island finder to define the island
	private final static boolean	PRINT_PROGRESS 		= false;	// set to true to print progress info
	private final static boolean	PRINT_FILES 		= true;		// set to true to print the bed files with the data from the simulation
	private final static int		GAUSSIAN_MV_WIDTH 	= 400000;	// moving window width of the gaussian smoothing
	private final static float		IF_MIN_WINDOW 		= 0.02f;	// island finder minimum window score parameter
	private final static int 		IF_GAP 				= 500;		// island finder gap parameter
	private final static int 		IF_MIN_LENGTH 		= 50;		// island finder island minimum length parameter

	private final File		outputDir;					// directory for the output data
	private final int 		islandSize;					// size of the islands to use in the simulation
	private final double 	percentageReadToAdd;		// percentage of reads to add in the S phase in the islands
	private final SCWList 	sList;						// s phase data
	private final SCWList 	g1List;						// g1 phase data


	/**
	 * Creates an instance of {@link SingleSimulation}
	 * @param outputDir directory for the output data
	 * @param islandSize size of the islands to use in the simulation
	 * @param percentageReadToAdd percentage of reads to add in the S phase in the islands
	 * @param sList s phase data
	 * @param g1List g1 phase data
	 */
	public SingleSimulation(File outputDir,
			int islandSize,
			double percentageReadToAdd,
			SCWList sList,
			SCWList g1List) {
		this.outputDir = outputDir;
		this.islandSize = islandSize;
		this.percentageReadToAdd = percentageReadToAdd;
		this.sList = sList;
		this.g1List = g1List;
	}


	@Override
	public SimulationResult compute() throws Exception {
		printProgress("SingleSimulation.compute() - 1");
		// 1 - generate control lists
		SCWList[] resampledList = new ResampleLayers(sList, g1List, 0).compute();
		SCWList controlS = resampledList[0];
		SCWList controlG1 = resampledList[1];

		// 2 - generate sample lists

		// 2a - generate list with reads added
		printProgress("SingleSimulation.compute() - 2a");
		resampledList = new ResampleLayers(sList, g1List, percentageReadToAdd).compute();
		SCWList resampledSReadAdded = resampledList[0];
		SCWList resampledG1ReadAdded = resampledList[1];

		// 2b - generate list with no reads added
		printProgress("SingleSimulation.compute() - 2b");
		resampledList = new ResampleLayers(sList, g1List, 0).compute();
		SCWList resampledSNoReadAdded = resampledList[0];
		SCWList resampledG1NoReadAdded = resampledList[1];

		// 2c - generate islands mask list and its inverse
		printProgress("SingleSimulation.compute() - 2c");
		SCWList islandMask = new GenerateIslands(ISLAND_DISTANCE, islandSize, g1List).compute();
		SCWList islandInvertedMask = new MCWLOInvertMask(islandMask).compute();
		//printSCWInTmpFile(islandMask.get(0), "inislands_" + islandSize + "_" + percentageReadToAdd + "_");

		// 2d - sum the island with reads added with the baseline
		printProgress("SingleSimulation.compute() - 2d");
		resampledSReadAdded = new SCWLOTwoLayers(resampledSReadAdded, islandMask, ScoreOperation.MULTIPLICATION).compute();
		resampledSNoReadAdded = new SCWLOTwoLayers(resampledSNoReadAdded, islandInvertedMask, ScoreOperation.MULTIPLICATION).compute();
		SCWList resampledS = new SCWLOTwoLayers(resampledSReadAdded, resampledSNoReadAdded, ScoreOperation.ADDITION).compute();

		resampledG1ReadAdded = new SCWLOTwoLayers(resampledG1ReadAdded, islandMask, ScoreOperation.MULTIPLICATION).compute();
		resampledG1NoReadAdded = new SCWLOTwoLayers(resampledG1NoReadAdded, islandInvertedMask, ScoreOperation.MULTIPLICATION).compute();
		SCWList resampledG1 = new SCWLOTwoLayers(resampledG1ReadAdded, resampledG1NoReadAdded, ScoreOperation.ADDITION).compute();

		// 3 - convert into binlist
		printProgress("SingleSimulation.compute() - 3");
		BinList binnedControlS = new SCWLOConvertIntoBinList(controlS, 500, ScoreOperation.ADDITION).compute();
		BinList binnedControlG1 = new SCWLOConvertIntoBinList(controlG1, 500, ScoreOperation.ADDITION).compute();
		BinList binnedResampledS = new SCWLOConvertIntoBinList(resampledS, 500, ScoreOperation.ADDITION).compute();
		BinList binnedResampledG1 = new SCWLOConvertIntoBinList(resampledG1, 500, ScoreOperation.ADDITION).compute();

		// 4 - gauss binlists
		printProgress("SingleSimulation.compute() - 4");
		binnedControlS = new BLOGauss(binnedControlS, GAUSSIAN_MV_WIDTH, false).compute();
		binnedControlG1 = new BLOGauss(binnedControlG1, GAUSSIAN_MV_WIDTH, false).compute();
		binnedResampledS = new BLOGauss(binnedResampledS, GAUSSIAN_MV_WIDTH, false).compute();
		binnedResampledG1 = new BLOGauss(binnedResampledG1, GAUSSIAN_MV_WIDTH, false).compute();

		// 5 - compute S / G1 ratios
		printProgress("SingleSimulation.compute() - 5");
		BinList controlSG1 = (BinList) new BLOTwoLayers(binnedControlS, binnedControlG1, ScoreOperation.DIVISION).compute();
		BinList sampleSG1 = (BinList) new BLOTwoLayers(binnedResampledS, binnedResampledG1, ScoreOperation.DIVISION).compute();

		// 6 - remove windows that are null in one of the 2 lists
		printProgress("SingleSimulation.compute() - 6");
		BinList maskControl = (BinList) new SCWLOOperationWithConstant(controlSG1, OperationWithConstant.UNIQUE_SCORE, 1, false).compute();
		BinList maskSample = (BinList) new SCWLOOperationWithConstant(sampleSG1, OperationWithConstant.UNIQUE_SCORE, 1, false).compute();
		controlSG1 = (BinList) new BLOTwoLayers(maskSample, controlSG1, ScoreOperation.MULTIPLICATION).compute();
		sampleSG1 = (BinList) new BLOTwoLayers(maskControl, sampleSG1, ScoreOperation.MULTIPLICATION).compute();
		printSCWInTmpFile(controlSG1.get(0), "controlSG1" );
		printSCWInTmpFile(sampleSG1.get(0), "sampleSG1");

		// 7 - compute sample - control difference
		printProgress("SingleSimulation.compute() - 7");
		BinList sampleCtrlDifference = (BinList) new BLOTwoLayers(sampleSG1, controlSG1, ScoreOperation.SUBTRACTION).compute();
		printSCWInTmpFile(sampleCtrlDifference.get(0), "difference");

		// 8 - call islands
		printProgress("SingleSimulation.compute() - 8");
		GeneList islands;
		if (USE_ISLAND_FINDER) {
			BinList positiveSampleCtrlDifference = (BinList) new SCWLOFilterThreshold(sampleCtrlDifference, 0, Float.POSITIVE_INFINITY, false).compute();
			BinList negativesSampleCtrlDifference = (BinList) new SCWLOOperationWithConstant(sampleCtrlDifference, OperationWithConstant.MULTIPLICATION, -1f, false).compute();
			negativesSampleCtrlDifference = (BinList) new SCWLOFilterThreshold(negativesSampleCtrlDifference, 0, Float.POSITIVE_INFINITY, false).compute();
			GeneList positiveIslands = findIslandUsingIslandFinder(positiveSampleCtrlDifference);
			GeneList negativeIslands = findIslandUsingIslandFinder(negativesSampleCtrlDifference);
			islands = new GLOMergeGeneLists(positiveIslands, negativeIslands).compute();
			islands = flattenGeneList(islands);
		} else {
			islands = new FindIslands(sampleCtrlDifference).compute();
		}

		// 9 - score islands
		printProgress("SingleSimulation.compute() - 9");
		GeneList controlIslandsS = new GLOScoreFromSCWList(islands, controlS, GeneScoreType.BASE_COVERAGE_SUM).compute();
		GeneList controlIslandsG1 = new GLOScoreFromSCWList(islands, controlG1, GeneScoreType.BASE_COVERAGE_SUM).compute();
		GeneList sampleIslandsS = new GLOScoreFromSCWList(islands, resampledS, GeneScoreType.BASE_COVERAGE_SUM).compute();
		GeneList sampleIslandsG1 = new GLOScoreFromSCWList(islands, resampledG1, GeneScoreType.BASE_COVERAGE_SUM).compute();

		// 10 - compute fisher exact test and retrieve qvalues
		printProgress("SingleSimulation.compute() - 10");
		SCWList islandsQValues = new ComputeQValues(controlIslandsS, controlIslandsG1, sampleIslandsS, sampleIslandsG1).compute();
		printSCWInTmpFile(islandsQValues.get(0), "islands");

		// 11 - filter islands with qvalue under 0.05
		printProgress("SingleSimulation.compute() - 11");
		SCWList filteredIslands = new SCWLOFilterThreshold(islandsQValues, 0, Q_VALUE_CUTOFF, false).compute();
		filteredIslands = new SCWLOConvertIntoSimpleSCWList(filteredIslands, SCWListType.MASK).compute();

		// 12 - compute false positive and false negatives
		printProgress("SingleSimulation.compute() - 12");
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
		bloFindIslands.getIsland().setWindowMinValue(IF_MIN_WINDOW);
		bloFindIslands.getIsland().setGap(IF_GAP);
		bloFindIslands.getIsland().setIslandMinScore(0);
		bloFindIslands.getIsland().setIslandMinLength(IF_MIN_LENGTH);
		IslandResultType[] resType = {IslandResultType.IFSCORE};
		bloFindIslands.setList(resType);
		BinList resBinList = bloFindIslands.compute()[0];
		SCWList resMaskList = new SCWLOConvertIntoSimpleSCWList(resBinList, SCWListType.MASK).compute();
		GeneList resGeneList = new SCWLOConvertIntoGeneList(resMaskList).compute();
		return resGeneList;
	}

	/**
	 * Merge overlapping genes of gene list together
	 * @param geneList a gene list
	 * @return a new gene list
	 * @throws Exception
	 */
	private GeneList flattenGeneList(GeneList geneList) throws Exception {
		SCWList mask = new SCWLOConvertIntoSimpleSCWList(geneList, SCWListType.MASK).compute();
		return new SCWLOConvertIntoGeneList(mask).compute();
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


	/**
	 * Prints the specified string in the standard output if {@link #PRINT_PROGRESS} is set to true
	 * @param stringToPrint
	 */
	private void printProgress(String stringToPrint) {
		if (PRINT_PROGRESS) {
			System.out.println(stringToPrint);
		}
	}


	/**
	 * Prints the specified {@link ListView} in a temporary file with the specified name prefix
	 * @param data
	 * @param prefix
	 * @throws IOException
	 */
	public void printSCWInTmpFile(ListView<? extends ScoredChromosomeWindow> data, String prefix) throws IOException {
		if (PRINT_FILES) {
			File file = new File(outputDir.getAbsolutePath() + islandSize + "_" + percentageReadToAdd + "_" + prefix + ".bed");
			System.out.println("Writing file: " + file.getPath());
			//Write the data
			PrintWriter dout = new PrintWriter(file);
			for (int i = 0; i < data.size(); i++) {
				ScoredChromosomeWindow scw = data.get(i);
				dout.println("chr1\t" + (scw.getStart() - 1)+ "\t" + (scw.getStop() - 1) + "\t-\t" + scw.getScore());
			}
			dout.close();
		}
	}


	@Override
	public void stop() {
		// stop operation not implemented
	}
}
