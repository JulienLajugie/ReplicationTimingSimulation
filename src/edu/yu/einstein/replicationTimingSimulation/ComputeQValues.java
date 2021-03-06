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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.security.InvalidParameterException;
import java.util.concurrent.ExecutionException;

import edu.yu.einstein.genplay.core.manager.project.ProjectChromosomes;
import edu.yu.einstein.genplay.core.manager.project.ProjectManager;
import edu.yu.einstein.genplay.core.operation.Operation;
import edu.yu.einstein.genplay.dataStructure.chromosome.Chromosome;
import edu.yu.einstein.genplay.dataStructure.list.chromosomeWideList.SCWListView.dense.DenseSCWListViewBuilder;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.SCWList.SCWList;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.SCWList.SCWListBuilder;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.geneList.GeneList;

/**
 * Computes the q-values from the p-values of the fisher exact test
 * of the ratios sample S / G1 and control S/ G1
 * @author Julien Lajugie
 */
public class ComputeQValues implements Operation<SCWList> {

	private final GeneList 	controlIslandsS;
	private final GeneList 	controlIslandsG1;
	private final GeneList 	sampleIslandsS;
	private final GeneList 	sampleIslandsG1;


	/**
	 * Creates an instance of {@link ComputeQValues}
	 * @param controlIslandsS
	 * @param controlIslandsG1
	 * @param sampleIslandsS
	 * @param sampleIslandsG1
	 */
	public ComputeQValues(GeneList controlIslandsS, GeneList controlIslandsG1,
			GeneList sampleIslandsS, GeneList sampleIslandsG1) {
		this.controlIslandsS = controlIslandsS;
		this.controlIslandsG1 = controlIslandsG1;
		this.sampleIslandsS = sampleIslandsS;
		this.sampleIslandsG1 = sampleIslandsG1;
	}


	@Override
	public SCWList compute() throws Exception {
		SCWList qValueScwList = performFisherExactTestAndRetrieveQValues();
		return qValueScwList;
	}


	@Override
	public String getDescription() {
		return "Operation: Compute Q-Values";
	}


	@Override
	public String getProcessingDescription() {
		return "Computing Q-Values";
	}


	@Override
	public int getStepCount() {
		return 1;
	}


	private SCWList performFisherExactTestAndRetrieveQValues() throws IOException, InterruptedException, CloneNotSupportedException, InvalidParameterException, ExecutionException {
		File tmpD = File.createTempFile("fish_in_", ".txt");
		File tmpR = File.createTempFile("fish_out_", ".txt");
		File tmpScript = File.createTempFile("fish", ".R");
		File tmpOut = File.createTempFile("fish_rout", ".txt");

		//data = matrix with columns a, b, c, d
		//test will be performed on each row
		PrintWriter sout = new PrintWriter(tmpScript);
		//sout.println("library(qvalue)");
		sout.println("d = read.table('" + tmpD.getAbsolutePath() + "')");
		sout.println("f = apply(d, 1, function(x) {chisq.test(rbind(c(x[1],x[2]), c(x[3],x[4])))})");
		sout.println("p = as.numeric(lapply(f, function(x) { x$p.value }))");
		sout.println("p[p > 1] = 1");
		sout.println("q <- p.adjust(p, method=\"fdr\")");
		//sout.println("q = qvalue(p)");
		sout.println("write.table(cbind(p, q), file='" + tmpR.getAbsolutePath() + "', col.names=F, row.names=F)");
		sout.close();

		//Write the data
		PrintWriter dout = new PrintWriter(tmpD);
		for (int i = 0; i < controlIslandsS.size(); i++) {
			for (int j = 0; j < controlIslandsS.size(i); j++) {
				long a = (long) sampleIslandsS.get(i, j).getScore();
				long b = (long) sampleIslandsG1.get(i, j).getScore();
				long c = (long) controlIslandsS.get(i, j).getScore();
				long d = (long) controlIslandsG1.get(i, j).getScore();
				if ((a != 0) || (b != 0) || (c != 0) || (d != 0)) {
					dout.println(a + "\t" + b + "\t" + c + "\t" + d);
				}
			}
		}
		dout.close();
		//Run the R command
		String cmd = "R CMD BATCH " + tmpScript.getAbsolutePath() + " " + tmpOut.getAbsolutePath();
		System.out.println(cmd);
		Runtime.getRuntime().exec(cmd).waitFor();
		//Read the output
		BufferedReader in = new BufferedReader(new FileReader(tmpR));
		DenseSCWListViewBuilder prototypeBuilder = new DenseSCWListViewBuilder();
		SCWListBuilder resultListBuilder = new SCWListBuilder(prototypeBuilder);
		ProjectChromosomes projectChromosomes = ProjectManager.getInstance().getProjectChromosomes();
		for (int i = 0; i < controlIslandsS.size(); i++) {
			Chromosome currentChromo = projectChromosomes.get(i);
			for (int j = 0; j < controlIslandsS.size(i); j++) {
				long a = (long) sampleIslandsS.get(i, j).getScore();
				long b = (long) sampleIslandsG1.get(i, j).getScore();
				long c = (long) controlIslandsS.get(i, j).getScore();
				long d = (long) controlIslandsG1.get(i, j).getScore();
				if ((a != 0) || (b != 0) || (c != 0) || (d != 0)) {
					// pValues in first column, qValues in second column
					String[] cols = in.readLine().split(" ", 2);
					double qValue = Double.parseDouble(cols[1]);
					// don't want the really small values to be rounded down to 0
					if (qValue < Float.MIN_NORMAL) {
						qValue = Float.MIN_NORMAL;
					}
					int start = controlIslandsS.get(i, j).getStart();
					int stop = controlIslandsS.get(i, j).getStop();
					resultListBuilder.addElementToBuild(currentChromo, start, stop, (float) qValue);
				}
			}
		}
		in.close();
		return resultListBuilder.getSCWList();
	}


	@Override
	public void stop() {
		// this operation cannot be stop
	}
}
