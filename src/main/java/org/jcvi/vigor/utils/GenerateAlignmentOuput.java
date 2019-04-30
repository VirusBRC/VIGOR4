package org.jcvi.vigor.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.Vigor;
import org.jcvi.vigor.component.Model;
import org.springframework.stereotype.Service;

@Service
public class GenerateAlignmentOuput {

	private static Logger LOGGER = LogManager.getLogger(Vigor.class);

	// public void generateOutputFile ( GenerateVigorOutput.Outfiles outfiles,
	// List<Model> models ) {
	public void generateOutputFile(String outputDir, GenerateVigorOutput.Outfiles outfiles, List<Model> models) {

		Model modelg = models.get(0);
		String GeneIDg = modelg.getGeneID();

		Pattern p = Pattern.compile(GenerateVigorOutput.REGEX);
		Pattern p2 = Pattern.compile(GenerateVigorOutput.REGEX2);

		Matcher mg = p2.matcher(GeneIDg);
		GeneIDg = mg.replaceAll(GenerateVigorOutput.REPLACE_G);

		String[] Parts = GeneIDg.split("-");
		String Part0 = Parts[0];

		Matcher mgg = p.matcher(Part0);
		GeneIDg = mgg.replaceAll(GenerateVigorOutput.REPLACE);

		String FileNameBRCg = GeneIDg.replaceAll("/", "_");

		String FilePathBRCg = outputDir;
		String FileNameg = FilePathBRCg + "/" + FileNameBRCg + ".aln";

		try {
			FileWriter writerg = new FileWriter(FileNameg);
			BufferedWriter bwBRCg = new BufferedWriter(writerg);

			List<File> raw_files = models.stream().map(m -> m.getAlignment().getAlignmentEvidence().getRaw_alignment()).distinct()
					.collect(Collectors.toList());

			for (File raw_alignment : raw_files) {
				// printAlignment(outfiles.get(GenerateVigorOutput.Outfile.ALN),
				// raw_alignment);
				printAlignment(bwBRCg, outfiles.get(GenerateVigorOutput.Outfile.ALN), raw_alignment);

			}

			List<File> temp_directories = models.stream().map(m -> m.getAlignment().getAlignmentEvidence().getResults_directory())
					.distinct().collect(Collectors.toList());
			for (File temp_directory : temp_directories) {
				VigorUtils.deleteTempFiles(temp_directory.getAbsolutePath().toString());
			}

			bwBRCg.close();
		} catch (Exception e) {
			System.out.println("Error in opening bwBRCg for outputing  aln" + e);
		}
	}

	// public void printAlignment ( BufferedWriter bw, File inputFile ) {
	public void printAlignment(BufferedWriter bwBRCg, BufferedWriter bw, File inputFile) {

		try {
			FileInputStream fRead = new FileInputStream(inputFile);
			int c;
			while ((c = fRead.read()) != -1) {
				bw.write((char) c);
				bwBRCg.write((char) c);
			}
			fRead.close();
		} catch (Exception e) {
			LOGGER.warn("Error reading temporory alignment file", inputFile.getAbsolutePath());
		}
	}
}
