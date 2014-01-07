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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.InvalidParameterException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import edu.yu.einstein.genplay.core.IO.dataReader.SCWReader;
import edu.yu.einstein.genplay.core.IO.extractor.Extractor;
import edu.yu.einstein.genplay.core.IO.extractor.ExtractorFactory;
import edu.yu.einstein.genplay.core.IO.genomeListLoader.AssemblyListLoader;
import edu.yu.einstein.genplay.core.manager.project.ProjectManager;
import edu.yu.einstein.genplay.dataStructure.chromosome.Chromosome;
import edu.yu.einstein.genplay.dataStructure.enums.ScoreOperation;
import edu.yu.einstein.genplay.dataStructure.enums.ScorePrecision;
import edu.yu.einstein.genplay.dataStructure.genome.Assembly;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.SCWList.SCWList;
import edu.yu.einstein.genplay.dataStructure.list.genomeWideList.SCWList.SCWListFactory;
import edu.yu.einstein.genplay.dataStructure.list.primitiveList.PrimitiveList;

/**
 * Class with main method to start the simulation
 * @author Julien Lajugie
 */
public class RunSimulationBatch {

	/**
	 * Class for JCommander to parse the command line arguments
	 * @author Julien Lajugie
	 */
	public static class Args {
		@Parameter(names = "-s", description = "File with replication timming data for the S phase")
		private String sFile;

		@Parameter(names = "-g1", description = "File with replication timming data for the G1 phase")
		private String g1File;

		@Parameter(names = "-out", description = "Output directory for the results of the simulation")
		private String outDir;
	}
	private final static int FALSE_POSITIVES = 0;
	private final static int FALSE_NEGATIVES = 1;
	private final static int ISLAND_AVERAGE_SIZE = 2;
	private final static int SAMPLE_CTRL_AVERAGE_DIFFERENCE = 3;

	// Island sizes to consider for the simulation
	private final static int[] islandSizes = {125000, 250000, 500000, 1000000, 2000000};

	// Percentage of reads to add for the simulation
	private final static double[] pctReadToAdds = {0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5};

	//private final static int[] islandSizes = {500000};
	//private final static double[] pctReadToAdds = {0.15};

	/**
	 * Initializes genplay project manager
	 * @throws ParserConfigurationException
	 * @throws SAXException
	 * @throws IOException
	 */
	public static void initManagers() throws ParserConfigurationException, SAXException, IOException {
		AssemblyListLoader listLoader = new AssemblyListLoader();
		// retrieve human assembly
		Assembly assembly = listLoader.getCladeList().get("mammal").getGenomeList().get("human").getAssemblyList().get("2009 02 hg19");
		// keep only the "basic" chromosomes (we don't want the _random)
		List<Chromosome> chrList = assembly.getChromosomeList();
		Collections.sort(chrList);
		chrList = chrList.subList(0, 23);
		assembly.setChromosomeList(chrList);
		ProjectManager.getInstance().setAssembly(assembly);
		ProjectManager.getInstance().updateChromosomeList();
		// set the precision of the data to 32 bit
		PrimitiveList.setScorePrecision(ScorePrecision.PRECISION_32BIT);
	}


	/**
	 * Extracts and generates a {@link SCWList} from the specified file
	 * @param file
	 * @return
	 * @throws Exception
	 */
	public static SCWList loadInputFile(File file) throws Exception {
		Extractor extractor = ExtractorFactory.getExtractor(file);
		SCWList scwList = SCWListFactory.createDenseSCWList((SCWReader) extractor, ScoreOperation.ADDITION);
		return scwList;
	}


	/**
	 * Main method, starts the simulation
	 * @param args
	 */
	public static void main(String[] args) {
		Args parameters = new Args();
		new JCommander(parameters, args);
		try {
			File sFile = new File(parameters.sFile);
			File g1File = new File(parameters.g1File);
			File outDir = new File(parameters.outDir);
			if (!outDir.exists()) {
				outDir.mkdir();
			}
			File outFile = new File(outDir, "simulation_summary.tsv");

			initManagers();

			SCWList sList = loadInputFile(sFile);
			SCWList g1List = loadInputFile(g1File);

			List<SimulationResult> resultList = new ArrayList<SimulationResult>();
			for (double pctReadToAdd: pctReadToAdds) {
				for (int islandSize: islandSizes) {
					System.out.println("*** Simulation with "
							+ (pctReadToAdd * 100)
							+ "% reads added on islands of "
							+ NumberFormat.getIntegerInstance().format(islandSize)
							+ "bp starting ***");
					SimulationResult result = new SingleSimulation(outDir, islandSize, pctReadToAdd, sList, g1List).compute();
					resultList.add(result);
				}
			}

			printResult(outFile, resultList);
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			System.exit(0);
		}
	}


	private static void printResult(BufferedWriter writer, List<SimulationResult> resultList, String title, int elementToPrint) throws IOException {
		writer.write(title);
		writer.newLine();
		for (int islandSize: islandSizes) {
			writer.write("\t" + islandSize);
		}
		writer.newLine();
		int index = 0;
		for (int i = 0; i < pctReadToAdds.length; i++) {
			writer.write(Double.toString(pctReadToAdds[i]));
			for (int j = 0; j < islandSizes.length; j++) {
				switch (elementToPrint) {
				case FALSE_POSITIVES:
					writer.write("\t" + resultList.get(index).getFalsePositiveRate());
					break;
				case FALSE_NEGATIVES:
					writer.write("\t" + resultList.get(index).getFalseNegativeRate());
					break;
				case ISLAND_AVERAGE_SIZE:
					writer.write("\t" + resultList.get(index).getIslandAverageSize());
					break;
				case SAMPLE_CTRL_AVERAGE_DIFFERENCE:
					writer.write("\t" + resultList.get(index).getSampleCtrlAverageDifference());
					break;
				default:
					throw new InvalidParameterException("Invalid element to print");
				}
				index++;
			}
			writer.newLine();
		}
		writer.newLine();
	}


	/**
	 * Prints the result of the simulation in the specified file
	 * @param outFile
	 * @param resultList
	 * @throws IOException
	 */
	private static void printResult(File outFile, List<SimulationResult> resultList) throws IOException {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(outFile));
			printResult(writer, resultList, "FALSE POSITIVES", FALSE_POSITIVES);
			printResult(writer, resultList, "FALSE NEGATIVES", FALSE_NEGATIVES);
			printResult(writer, resultList, "ISLAND AVERAGE SIZE", ISLAND_AVERAGE_SIZE);
			printResult(writer, resultList, "SAMPLE - CTRL AVERAGE DIFFERENCE", SAMPLE_CTRL_AVERAGE_DIFFERENCE);
		} finally {
			if (writer != null) {
				writer.close();
			}
		}
	}
}
