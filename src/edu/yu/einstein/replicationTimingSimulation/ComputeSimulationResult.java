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

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;

import edu.yu.einstein.genplay.core.manager.project.ProjectChromosomes;
import edu.yu.einstein.genplay.core.manager.project.ProjectManager;
import edu.yu.einstein.genplay.core.operation.Operation;
import edu.yu.einstein.genplay.core.operationPool.OperationPool;
import edu.yu.einstein.genplay.dataStructure.chromosome.Chromosome;
import edu.yu.einstein.genplay.dataStructure.chromosomeWindow.ChromosomeWindow;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.SCWList.SCWList;
import edu.yu.einstein.genplay.dataStructure.list.listView.ListView;
import edu.yu.einstein.genplay.dataStructure.scoredChromosomeWindow.ScoredChromosomeWindow;
import edu.yu.einstein.genplay.util.ListView.ListViews;

/**
 * Creates an instance of {@link ComputeSimulationResult}.
 * Compute the number of false positives and false negatives.
 * The false positives are the islands that are in the islands found but not in the islands generated
 * The false negatives are the islands that are in the islands generated but not in the islands found
 **/
public class ComputeSimulationResult implements Operation<SimulationResult> {

	private final int 		islandSize;						// size of the islands used in the simulation
	private final double 	percentageReadsAdded;			// number of reads added to the island (eg: 0.1 if there were 10% more reads)
	private final SCWList 	islandMasks;					// scw list containing the island generated
	private final SCWList 	islandsFound;					// scw list containing the island found during the simulation
	private final double	sampleCtrlAverageDifference;	// average difference between S and G1 on the island after gaussing
	private final double	sampleCtrlDifferenceStdErr;		// standard error difference between the sample on the islands after gaussing
	private boolean			stopped = false;				// true if the operation must be stopped


	/**
	 * Creates an instance of {@link ComputeSimulationResult}.
	 * Compute the number of false positives and false negatives.
	 * The false positives are the islands that are in the islands found but not in the islands generated
	 * The false negatives are the islands that are in the islands generated but not in the islands found
	 * @param islandSize size of the islands used in the simulation
	 * @param percentageReadsAdded number of reads added to the island (eg: 0.1 if there were 10% more reads)
	 * @param islandMasks mask containing the island generated
	 * @param islandsFound gene list containing the island found during the simulation
	 * @param sampleCtrlAverageDifference average difference between the sample on the islands after gaussing
	 * @param sampleCtrlDifferenceStdErr standard error difference between the sample on the islands after gaussing
	 */
	public ComputeSimulationResult(int islandSize, double percentageReadsAdded,
			SCWList islandMasks, SCWList islandsFound,
			double sampleCtrlAverageDifference, double  sampleCtrlDifferenceStdErr) {
		this.islandSize = islandSize;
		this.percentageReadsAdded = percentageReadsAdded;
		this.islandMasks = islandMasks;
		this.islandsFound = islandsFound;
		this.sampleCtrlAverageDifference = sampleCtrlAverageDifference;
		this.sampleCtrlDifferenceStdErr = sampleCtrlDifferenceStdErr;
	}


	@Override
	public SimulationResult compute() throws Exception {
		int islandCreatedCount = (int) islandMasks.getStatistics().getWindowCount();
		int islandFoundCount = (int) islandsFound.getStatistics().getWindowCount();

		int falsePositiveCount = 0;
		int falseNegativeCount = 0;

		ProjectChromosomes projectChromosomes = ProjectManager.getInstance().getProjectChromosomes();
		final OperationPool op = OperationPool.getInstance();
		final Collection<Callable<int[]>> threadList = new ArrayList<Callable<int[]>>();

		for(final Chromosome currentChromosome : projectChromosomes) {
			final ListView<ScoredChromosomeWindow> currentGeneratedIslands = islandMasks.get(currentChromosome);
			final ListView<ScoredChromosomeWindow> currentFoundIslands = islandsFound.get(currentChromosome);
			Callable<int[]> currentThread = new Callable<int[]>() {

				@Override
				public int[] call() throws Exception {
					// array for the result, the first element corresponds the false positives, the second to the false negatives
					int[] result = new int[2];
					// count false positives
					for (ChromosomeWindow foundWindow: currentFoundIslands) {
						if (stopped) {
							return null;
						}
						if (!hasOverlap(foundWindow, currentGeneratedIslands)) {
							result[0]++;
						}
					}
					if (stopped) {
						return null;
					}
					// count false negatives
					for (ChromosomeWindow generatedWindow: currentGeneratedIslands) {
						if (stopped) {
							return null;
						}
						if (!hasOverlap(generatedWindow, currentFoundIslands)) {
							result[1]++;
						}
					}
					op.notifyDone();
					return result;
				}
			};
			threadList.add(currentThread);
		}
		List<int[]> result = op.startPool(threadList);
		// sum up false positives and false negatives
		for (int[] currentResult: result) {
			falsePositiveCount += currentResult[0];
			falseNegativeCount += currentResult[1];
		}
		// compute the average size of the islands
		double islandSizeSum = islandsFound.getStatistics().getWindowLength();
		int islandAverageSize = 0;
		if (islandFoundCount != 0) {
			islandAverageSize = (int) (islandSizeSum / islandFoundCount);
		}
		double islandSizeStdErr = computeIslandSizeStdErr(islandFoundCount, islandAverageSize);
		return new SimulationResult(islandSize, percentageReadsAdded, islandCreatedCount, islandFoundCount, falsePositiveCount, falseNegativeCount, islandAverageSize, islandSizeStdErr, sampleCtrlAverageDifference, sampleCtrlDifferenceStdErr);
	}


	/**
	 * @param islandFoundCount number of islands found in the simulation
	 * @param islandAverageSize average size of the islands  found in the simulation
	 * @return the standard error of the length of the islands found in the simulation
	 * @throws InterruptedException
	 * @throws ExecutionException
	 */
	private double computeIslandSizeStdErr(int islandFoundCount, final double islandAverageSize) throws InterruptedException, ExecutionException {
		// compute the std deviation of the island length
		ProjectChromosomes projectChromosomes = ProjectManager.getInstance().getProjectChromosomes();
		final OperationPool op = OperationPool.getInstance();
		final Collection<Callable<Double>> threadList = new ArrayList<Callable<Double>>();

		for(final Chromosome currentChromosome : projectChromosomes) {
			final ListView<ScoredChromosomeWindow> currentFoundIslands = islandsFound.get(currentChromosome);
			Callable<Double> currentThread = new Callable<Double>() {

				@Override
				public Double call() throws Exception {
					double result = 0;
					// compute standard error
					for (ChromosomeWindow generatedWindow: currentFoundIslands) {
						if (stopped) {
							return null;
						}
						result += Math.pow(generatedWindow.getSize() - islandAverageSize, 2);
						}
					op.notifyDone();
					return result;
				}
			};
			threadList.add(currentThread);
		}
		List<Double> result = op.startPool(threadList);
		double islandStdDev = 0;
		// sum up false positives and false negatives
		for (Double currentResult: result) {
			islandStdDev += currentResult;
		}
		islandStdDev = Math.sqrt(islandStdDev);
		// compute the standard error of the island length
		return islandStdDev / Math.sqrt(islandFoundCount);
	}


	@Override
	public String getDescription() {
		return "Operation: Compute Simulation Result";
	}


	@Override
	public String getProcessingDescription() {
		return "Computing Simulation Result";
	}


	@Override
	public int getStepCount() {
		return 1;
	}


	/**
	 * @param window
	 * @param windowList
	 * @return true if at least one element of the window list overlap the specified window. False otherwise
	 */
	private boolean hasOverlap(ChromosomeWindow window, ListView<? extends ChromosomeWindow> windowList) {
		if (windowList.isEmpty()) {
			return false;
		}
		int indexFound = ListViews.binarySearch(windowList, window);
		if (indexFound >= 0) {
			return true;
		} else {
			int indexBefore = -indexFound - 2;

			int indexAfter = -indexFound - 1;
			if (indexBefore < 0) {
				return windowList.get(indexAfter).getStart() < window.getStop();
			} else if (indexAfter >= windowList.size()) {
				return windowList.get(indexBefore).getStop() > window.getStart();
			} else {
				return (windowList.get(indexBefore).getStop() > window.getStart()) ||
						(windowList.get(indexAfter).getStart() < window.getStop());
			}
		}
	}


	@Override
	public void stop() {
		stopped = true;
	}
}
