package edu.yu.einstein.replicationTimingSimulation;

import edu.yu.einstein.genplay.core.operation.Operation;
import edu.yu.einstein.genplay.core.operation.SCWList.SCWLOConvertIntoBinList;
import edu.yu.einstein.genplay.core.operation.binList.BLOIntervalsScoring;
import edu.yu.einstein.genplay.dataStructure.enums.ScoreOperation;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.SCWList.SCWList;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.SCWList.binList.BinList;


/**
 * Computes the average value of the differences between the sample and the control on the islands
 * @author Julien Lajugie
 */
public class ComputeSampleCtrlAverageDifference implements Operation<Float>{

	private final BinList 	differenceList;		// list containing the differences between the sample and the control
	private final SCWList 	filteredIslands;	// islands found during the simulation
	private boolean			stopped = false;	// true if the operation must be stopped


	/**
	 * Creates an instance of {@link ComputeSampleCtrlAverageDifference}
	 * @param differenceList list containing the differences between the sample and the control
	 * @param islands islands found during the simulation
	 */
	public ComputeSampleCtrlAverageDifference(BinList differenceList, SCWList filteredIslands) {
		this.differenceList = differenceList;
		this.filteredIslands = filteredIslands;
	}


	@Override
	public Float compute() throws Exception {
		BinList filteredIslandsBinList = new SCWLOConvertIntoBinList(filteredIslands, 500, ScoreOperation.ADDITION).compute();
		if (stopped) {
			return null;
		}
		BinList averageList = new BLOIntervalsScoring(filteredIslandsBinList, differenceList, 100, ScoreOperation.AVERAGE).compute();
		return (float) averageList.getStatistics().getAverage();
	}


	@Override
	public String getDescription() {
		return "Operation: Compute Sample - Control Average Difference";
	}


	@Override
	public String getProcessingDescription() {
		return "Computing the Average Difference Between the Sample and the Control";
	}


	@Override
	public int getStepCount() {
		return 2;
	}


	@Override
	public void stop() {
		stopped = true;
	}
}
