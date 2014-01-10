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

/**
 * Represents the result of a simulation
 * @author Julien Lajugie
 */
public class SimulationResult {

	private final int 		islandSize;						// size of the islands used in the simulation
	private final double 	percentageReadsAdded;			// number of reads added to the island (eg: 0.1 if there were 10% more reads)
	private final int 		islandCreatedCount;				// number of islands created for the simulation
	private final int 		islandFoundCount;				// number of island found at the end of the simulation
	private final int 		falsePositiveCount;				// number of false positives in the simulation
	private final int 		falseNegativeCount;				// number of false negative in the simulation
	private final int 		islandAverageSize;				// average size of the island found
	private final double	islandSizeStdErr;				// standard error of the island size
	private final double	sampleCtrlAverageDifference;	// average difference between the sample and the control after gaussing
	private final double	sampleCtrlDifferenceStdErr;		// standard error of the differences between the sample and the control after gaussing


	/**
	 * Creates an instance of {@link SimulationResult}
	 * @param islandSize size of the islands used in the simulation
	 * @param percentageReadsAdded number of reads added to the island (eg: 0.1 if there were 10% more reads)
	 * @param islandCreatedCount number of island generated for the simulation
	 * @param islandFoundCount number of island detected during the simulation
	 * @param falsePositiveCount number of islands found that were not generated
	 * @param falseNegativeCount number of islands missed
	 * @param islandAverageSize average size of the island found during the simulation
	 * @param SG1AverageDifference average difference between the sample and the control after gaussing
	 */
	public SimulationResult(int islandSize, double percentageReadsAdded,
			int islandCreatedCount, int islandFoundCount,
			int falsePositiveCount, int falseNegativeCount,
			int islandAverageSize, double islandSizeStdErr,
			double sampleCtrlAverageDifference, double sampleCtrlDifferenceStdErr) {
		this.islandSize = islandSize;
		this.percentageReadsAdded = percentageReadsAdded;
		this.islandCreatedCount = islandCreatedCount;
		this.islandFoundCount = islandFoundCount;
		this.falsePositiveCount = falsePositiveCount;
		this.falseNegativeCount = falseNegativeCount;
		this.islandAverageSize = islandAverageSize;
		this.islandSizeStdErr = islandSizeStdErr;
		this.sampleCtrlAverageDifference = sampleCtrlAverageDifference;
		this.sampleCtrlDifferenceStdErr = sampleCtrlDifferenceStdErr;
	}


	/**
	 * @return the number of false negatives
	 */
	public int getFalseNegativeCount() {
		return falseNegativeCount;
	}


	/**
	 * @return the false negative rate of the simulation
	 */
	public double getFalseNegativeRate() {
		return falseNegativeCount / (double) islandCreatedCount;
	}


	/**
	 * @return the number of false positives
	 */
	public int getFalsePositiveCount() {
		return falsePositiveCount;
	}


	/**
	 * @return the false positive rate of the simulation
	 */
	public double getFalsePositiveRate() {
		if (islandFoundCount == 0) {
			return 0;
		} else {
			return falsePositiveCount / (double) islandFoundCount;
		}
	}


	/**
	 * @return the average size of the island found during the simulation
	 */
	public int getIslandAverageSize() {
		return islandAverageSize;
	}


	/**
	 * @return the number of islands created for the simulation
	 */
	public int getIslandCreatedCount() {
		return islandCreatedCount;
	}


	/**
	 * @return the number of island found at the end of the simulation
	 */
	public int getIslandFoundCount() {
		return islandFoundCount;
	}


	/**
	 * @return the size of the island used during the simulation
	 */
	public int getIslandSize() {
		return islandSize;
	}


	/**
	 * @return the standard error of the island size
	 */
	public double getIslandSizeStdErr() {
		return islandSizeStdErr;
	}


	/**
	 * @return the percentage of read added to the S phase in the island during the simulation
	 */
	public double getPercentageReadsAdded() {
		return percentageReadsAdded;
	}


	/**
	 * @return the average difference (after gaussing) between the sample and the control after gaussing
	 */
	public double getSampleCtrlAverageDifference() {
		return sampleCtrlAverageDifference;
	}


	/**
	 * @return the standard error of the differences between the sample and the control
	 */
	public double getSampleCtrlDifferenceStdErr() {
		return sampleCtrlDifferenceStdErr;
	}
}
