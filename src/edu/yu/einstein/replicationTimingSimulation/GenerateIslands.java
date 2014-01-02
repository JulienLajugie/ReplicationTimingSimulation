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

import edu.yu.einstein.genplay.core.manager.project.ProjectChromosomes;
import edu.yu.einstein.genplay.core.manager.project.ProjectManager;
import edu.yu.einstein.genplay.core.operation.Operation;
import edu.yu.einstein.genplay.core.operation.SCWList.SCWLOConvertIntoSimpleSCWList;
import edu.yu.einstein.genplay.core.operation.geneList.GLOFilterThreshold;
import edu.yu.einstein.genplay.core.operation.geneList.GLOScoreFromSCWList;
import edu.yu.einstein.genplay.core.operationPool.OperationPool;
import edu.yu.einstein.genplay.dataStructure.chromosome.Chromosome;
import edu.yu.einstein.genplay.dataStructure.enums.GeneScoreType;
import edu.yu.einstein.genplay.dataStructure.enums.SCWListType;
import edu.yu.einstein.genplay.dataStructure.enums.Strand;
import edu.yu.einstein.genplay.dataStructure.gene.Gene;
import edu.yu.einstein.genplay.dataStructure.list.chromosomeWideList.geneListView.GeneListViewBuilder;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.SCWList.SCWList;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.SCWList.SimpleSCWList.SimpleSCWList;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.geneList.GeneList;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.geneList.SimpleGeneList;
import edu.yu.einstein.genplay.dataStructure.list.listView.ListView;
import edu.yu.einstein.genplay.util.ListView.SCWListViews;

/**
 * Generates the island used during the simulation
 * @author Julien Lajugie
 */
public class GenerateIslands implements Operation<SCWList> {

	private final int 		stepSize;		// size between two islands start positions
	private final int 		islandSize;		// size of the islands
	private final SCWList	mappableData;	// input data, islands not overlapping data are removed (correspond to unmappable regions of the genome)
	private boolean			stopped = false;// true if the operation must be stopped


	/**
	 * Creates an instance of {@link GenerateIslands}
	 * @param stepSize size between two islands start positions
	 * @param islandSize size of the islands
	 * @param mappleData input data, islands not overlapping data are removed (correspond to unmappable regions of the genome)
	 */
	public GenerateIslands(int stepSize, int islandSize, SCWList mappleData) {
		this.stepSize = stepSize;
		this.islandSize = islandSize;
		mappableData = mappleData;
	}


	@Override
	public SCWList compute() throws Exception {
		final OperationPool op = OperationPool.getInstance();
		final Collection<Callable<ListView<Gene>>> threadList = new ArrayList<Callable<ListView<Gene>>>();
		ProjectChromosomes projectChromosomes = ProjectManager.getInstance().getProjectChromosomes();

		for (final Chromosome currentChr: projectChromosomes) {
			Callable<ListView<Gene>> currentThread = new Callable<ListView<Gene>>() {
				@Override
				public ListView<Gene> call() throws Exception {
					GeneListViewBuilder lvBuilder = new GeneListViewBuilder();
					int islandStart = 1;
					int islandStop = islandStart + islandSize;
					while ((islandStop < currentChr.getLength()) && !stopped) {
						lvBuilder.addElementToBuild("", Strand.FIVE, islandStart, islandStop, 1, islandStart, islandStop, SCWListViews.createGenericSCWListView(islandStart, islandStop, 1));
						islandStart += stepSize;
						islandStop = islandStart + islandSize;
					}
					// tell the operation pool that a chromosome is done
					op.notifyDone();
					return lvBuilder.getListView();
				}
			};

			threadList.add(currentThread);
		}
		List<ListView<Gene>> result = op.startPool(threadList);
		if (result != null) {
			GeneList resultList = new SimpleGeneList(result, null, null);
			// score the exon of the list
			resultList = new GLOScoreFromSCWList(resultList, mappableData, GeneScoreType.MAXIMUM_COVERAGE).compute();
			// remove genes where the score is 0
			resultList = new GLOFilterThreshold(resultList, Float.MIN_NORMAL, Float.POSITIVE_INFINITY, false).compute();
			// convert in mask list
			SCWList resultMaskList = new SCWLOConvertIntoSimpleSCWList(resultList, SCWListType.MASK).compute();
			return resultMaskList;
		} else {
			return null;
		}
	}


	@Override
	public String getDescription() {
		return "Operation: Generate Islands";
	}


	@Override
	public String getProcessingDescription() {
		return "Generating Islands";
	}


	@Override
	public int getStepCount() {
		return 3 + SimpleSCWList.getCreationStepCount(SCWListType.MASK);
	}


	@Override
	public void stop() {
		stopped = true;
	}
}
