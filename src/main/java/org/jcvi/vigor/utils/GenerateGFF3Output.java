package org.jcvi.vigor.utils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.VirusGenome;
import org.springframework.stereotype.Service;

@Service
public class GenerateGFF3Output {

	private static Range.CoordinateSystem oneBased = Range.CoordinateSystem.RESIDUE_BASED;

	// public void generateOutputFile ( GenerateVigorOutput.Outfiles outfiles,
	// List<Model> models ) throws IOException {
	public void generateOutputFile(String outputDir, GenerateVigorOutput.Outfiles outfiles, List<Model> models)
			throws IOException {

		// printGFF3Features(outfiles.get(GenerateVigorOutput.Outfile.GFF3),
		// models);
		printGFF3Features(outputDir, outfiles.get(GenerateVigorOutput.Outfile.GFF3), models);
	}

	// public void printGFF3Features ( BufferedWriter bw, List<Model> geneModels
	// ) throws IOException {
	public void printGFF3Features(String outputDir, BufferedWriter bw, List<Model> geneModels) throws IOException {

		Model modelg = geneModels.get(0);
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
		String FileNameg = FilePathBRCg + "/" + FileNameBRCg + ".gff3";

		FileWriter writerg = new FileWriter(FileNameg);
		BufferedWriter bwBRCg = new BufferedWriter(writerg);

		// Note: this string should be consistent with the one used in
		// Vigor.java
		bwBRCg.write("##gff-version 3\n");

		String idG;

		for (Model geneModel : geneModels) {
			List<String> notes = geneModel.getNotes();
			int i = 1;
			VirusGenome virusGenome = geneModel.getAlignment().getVirusGenome();
			long seqlength = virusGenome.getSequence().getLength();
			String geneomeSeqID = virusGenome.getId();
			List<Exon> exons = geneModel.getExons();
			String geneName = geneModel.getAlignment().getViralProtein().getGeneSymbol();
			long CDSStart = VigorFunctionalUtils.getDirectionBasedCoordinate(exons.get(0).getRange().getBegin(oneBased),
					seqlength, geneModel.getDirection());
			long CDSEnd = VigorFunctionalUtils.getDirectionBasedCoordinate(
					exons.get(exons.size() - 1).getRange().getEnd(oneBased), seqlength, geneModel.getDirection());
			String mRnaID = geneModel.getGeneID() + "." + i;

			bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
			bwBRCg.write(geneomeSeqID + "\t" + "vigor" + "\t");
			// gene
			if (geneModel.isPseudogene()) {
				bw.write("pseudogene" + "\t");
				bwBRCg.write("pseudogene" + "\t");
			} else {
				bw.write("gene" + "\t");
				bwBRCg.write("gene" + "\t");
			}
			bw.write(Math.min(CDSStart, CDSEnd) + "\t" + Math.max(CDSStart, CDSEnd) + "\t");
			bwBRCg.write(Math.min(CDSStart, CDSEnd) + "\t" + Math.max(CDSStart, CDSEnd) + "\t");
			bw.write("." + "\t");
			bwBRCg.write("." + "\t");
			if (geneModel.getDirection().equals(Direction.FORWARD)) {
				bw.write("+" + "\t");
				bwBRCg.write("+" + "\t");
			} else {
				bw.write("-" + "\t");
				bwBRCg.write("-" + "\t");
			}
			bw.write(exons.get(0).getFrame().getFrame() - 1 + "\t");
			bwBRCg.write(exons.get(0).getFrame().getFrame() - 1 + "\t");
			bw.write("ID=" + geneModel.getGeneID() + ";" + "Name=" + geneName + ";");
			bwBRCg.write("ID=" + geneModel.getGeneID() + ";" + "Name=" + geneName + ";");
			if (geneModel.isPartial3p() || geneModel.isPartial5p()) {
				bw.write("Partial" + ",");
				bwBRCg.write("Partial" + ",");
			}
			if (notes.contains(NoteType.Gene.toString())) {
				bw.write(String.format("Note=%s;", NoteType.Gene));
				bwBRCg.write(String.format("Note=%s;", NoteType.Gene));
			}
			bw.write("\n");
			bwBRCg.write("\n");
			// mRNA
			bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
			bwBRCg.write(geneomeSeqID + "\t" + "vigor" + "\t");
			bw.write("mRNA" + "\t");
			bwBRCg.write("mRNA" + "\t");
			bw.write(Math.min(CDSStart, CDSEnd) + "\t" + Math.max(CDSStart, CDSEnd) + "\t");
			bwBRCg.write(Math.min(CDSStart, CDSEnd) + "\t" + Math.max(CDSStart, CDSEnd) + "\t");
			bw.write("." + "\t");
			bwBRCg.write("." + "\t");
			if (geneModel.getDirection().equals(Direction.FORWARD)) {
				bw.write("+" + "\t");
				bwBRCg.write("+" + "\t");
			} else {
				bw.write("-" + "\t");
				bwBRCg.write("-" + "\t");
			}
			bw.write(exons.get(0).getFrame().getFrame() - 1 + "\t");
			bwBRCg.write(exons.get(0).getFrame().getFrame() - 1 + "\t");
			bw.write("ID=" + mRnaID + ";" + "Parent=" + geneModel.getGeneID() + ";");
			bwBRCg.write("ID=" + mRnaID + ";" + "Parent=" + geneModel.getGeneID() + ";");
			if (geneModel.isPartial3p() || geneModel.isPartial5p()) {
				bw.write("Partial;");
				bwBRCg.write("Partial;");
			}
			bw.write("\n");
			bwBRCg.write("\n");
			// exon
			IDGenerator idGenerator = new IDGenerator(mRnaID);
			for (int j = 0; j < exons.size(); j++) {
				Exon exon = exons.get(j);
				bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
				bwBRCg.write(geneomeSeqID + "\t" + "vigor" + "\t");
				bw.write("exon" + "\t");
				bwBRCg.write("exon" + "\t");
				long exonBegin = VigorFunctionalUtils.getDirectionBasedCoordinate(exon.getRange().getBegin(oneBased), seqlength,
						geneModel.getDirection());
				long exonEnd = VigorFunctionalUtils.getDirectionBasedCoordinate(exon.getRange().getEnd(oneBased), seqlength,
						geneModel.getDirection());
				bw.write(Math.min(exonBegin, exonEnd) + "\t" + Math.max(exonBegin, exonEnd) + "\t");
				bwBRCg.write(Math.min(exonBegin, exonEnd) + "\t" + Math.max(exonBegin, exonEnd) + "\t");
				bw.write("." + "\t");
				bwBRCg.write("." + "\t");
				bw.write((geneModel.getDirection() == Direction.FORWARD ? "+" : '-') + "\t");
				bwBRCg.write((geneModel.getDirection() == Direction.FORWARD ? "+" : '-') + "\t");
				bw.write(exon.getFrame().getFrame() - 1 + "\t");
				bwBRCg.write(exon.getFrame().getFrame() - 1 + "\t");
				idG = idGenerator.next();
				// bw.write(String.format("ID=%s;Parent=%s;",
				// idGenerator.next(), mRnaID));
				bw.write(String.format("ID=%s;Parent=%s;", idG, mRnaID));
				bwBRCg.write(String.format("ID=%s;Parent=%s;", idG, mRnaID));
				if ((j == 0 && geneModel.isPartial5p()) || (j == exons.size() - 1 && geneModel.isPartial3p())) {
					bw.write("Partial" + ";");
					bwBRCg.write("Partial" + ";");
				}
				bw.write("\n");
				bwBRCg.write("\n");
			}
			// CDS
			bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
			bwBRCg.write(geneomeSeqID + "\t" + "vigor" + "\t");
			bw.write("CDS" + "\t");
			bwBRCg.write("CDS" + "\t");
			bw.write(Math.min(CDSStart, CDSEnd) + "\t" + Math.max(CDSStart, CDSEnd) + "\t");
			bwBRCg.write(Math.min(CDSStart, CDSEnd) + "\t" + Math.max(CDSStart, CDSEnd) + "\t");
			bw.write("." + "\t");
			bwBRCg.write("." + "\t");
			if (geneModel.getDirection().equals(Direction.FORWARD)) {
				bw.write("+" + "\t");
				bwBRCg.write("+" + "\t");
			} else {
				bw.write("-" + "\t");
				bwBRCg.write("-" + "\t");
			}
			bw.write(exons.get(0).getFrame().getFrame() - 1 + "\t");
			bwBRCg.write(exons.get(0).getFrame().getFrame() - 1 + "\t");

			idG = idGenerator.next();
			// bw.write(String.format("ID=%s;Parent=%s;", idGenerator.next(),
			// geneModel.getGeneID()));
			bw.write(String.format("ID=%s;Parent=%s;", idG, geneModel.getGeneID()));
			bwBRCg.write(String.format("ID=%s;Parent=%s;", idG, geneModel.getGeneID()));
			if (geneModel.isPartial3p() || geneModel.isPartial5p()) {
				bw.write("Partial" + ";");
				bwBRCg.write("Partial" + ";");
			}
			bw.write("\n");
			bwBRCg.write("\n");
			// insertion
			if (geneModel.getInsertRNAEditingRange() != null) {
				bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
				bwBRCg.write(geneomeSeqID + "\t" + "vigor" + "\t");
				bw.write("insertion" + "\t");
				bwBRCg.write("insertion" + "\t");
				long insertBegin = VigorFunctionalUtils.getDirectionBasedCoordinate(
						geneModel.getInsertRNAEditingRange().getBegin(oneBased), seqlength, geneModel.getDirection());
				long insertEnd = VigorFunctionalUtils.getDirectionBasedCoordinate(
						geneModel.getInsertRNAEditingRange().getEnd(oneBased), seqlength, geneModel.getDirection());
				bw.write(Math.min(insertBegin, insertEnd) + "\t" + Math.min(insertBegin, insertEnd) + "\t");
				bwBRCg.write(Math.min(insertBegin, insertEnd) + "\t" + Math.min(insertBegin, insertEnd) + "\t");
				bw.write("." + "\t");
				bwBRCg.write("." + "\t");
				if (geneModel.getDirection().equals(Direction.FORWARD)) {
					bw.write("+" + "\t");
					bwBRCg.write("+" + "\t");
				} else {
					bw.write("-" + "\t");
					bwBRCg.write("-" + "\t");
				}
				bw.write("." + "\t");
				bwBRCg.write("." + "\t");

				idG = idGenerator.next();
				// bw.write(String.format("ID=%s;Parent=%s;",
				// idGenerator.next(), mRnaID));
				bw.write(String.format("ID=%s;Parent=%s;", idG, mRnaID));
				bwBRCg.write(String.format("ID=%s;Parent=%s;", idG, mRnaID));
				if (notes.contains(NoteType.RNA_Editing.toString())) {
					bw.write(String.format("Note=%s;", NoteType.RNA_Editing));
					bwBRCg.write(String.format("Note=%s;", NoteType.RNA_Editing));
				}
				bw.write("\n");
				bwBRCg.write("\n");
			}
			// Stop_codon_read_through
			if (geneModel.getReplaceStopCodonRange() != null) {
				bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
				bwBRCg.write(geneomeSeqID + "\t" + "vigor" + "\t");
				bw.write("stop_codon_read_through" + "\t");
				bwBRCg.write("stop_codon_read_through" + "\t");
				long replaceStopBegin = VigorFunctionalUtils.getDirectionBasedCoordinate(
						geneModel.getReplaceStopCodonRange().getBegin(oneBased), seqlength, geneModel.getDirection());
				long replaceStopEnd = VigorFunctionalUtils.getDirectionBasedCoordinate(
						geneModel.getReplaceStopCodonRange().getEnd(oneBased), seqlength, geneModel.getDirection());
				bw.write(Math.min(replaceStopBegin, replaceStopEnd) + "\t" + Math.max(replaceStopBegin, replaceStopEnd) + "\t");
				bwBRCg.write(
						Math.min(replaceStopBegin, replaceStopEnd) + "\t" + Math.max(replaceStopBegin, replaceStopEnd) + "\t");
				bw.write("." + "\t");
				bwBRCg.write("." + "\t");
				if (geneModel.getDirection().equals(Direction.FORWARD)) {
					bw.write("+" + "\t");
					bwBRCg.write("+" + "\t");
				} else {
					bw.write("-" + "\t");
					bwBRCg.write("-" + "\t");
				}
				bw.write("." + "\t");
				bwBRCg.write("." + "\t");
				idG = idGenerator.next();
				// bw.write(String.format("ID=%s;Parent=%s;",
				// idGenerator.next(), mRnaID));
				bw.write(String.format("ID=%s;Parent=%s;", idG, mRnaID));
				bwBRCg.write(String.format("ID=%s;Parent=%s;", idG, mRnaID));
				if (notes.contains(NoteType.StopCodonReadThrough.toString())) {
					bw.write(String.format("Note=%s;", (NoteType.StopCodonReadThrough)));
					bwBRCg.write(String.format("Note=%s;", (NoteType.StopCodonReadThrough)));
				}
				bw.write("\n");
				bwBRCg.write("\n");
			}
			// plus/minus_1_translationally_frameshifted
			if (geneModel.getAlignment().getViralProtein().getGeneAttributes().getRibosomal_slippage().isHas_ribosomal_slippage()
					&& geneModel.getRibosomalSlippageRange() != null) {
				int frameshift = geneModel.getAlignment().getViralProtein().getGeneAttributes().getRibosomal_slippage()
						.getSlippage_frameshift();
				long slippageBegin = VigorFunctionalUtils.getDirectionBasedCoordinate(
						geneModel.getRibosomalSlippageRange().getBegin(oneBased), seqlength, geneModel.getDirection());
				long slippageEnd = VigorFunctionalUtils.getDirectionBasedCoordinate(
						geneModel.getRibosomalSlippageRange().getEnd(oneBased), seqlength, geneModel.getDirection());
				if (frameshift == -1 || frameshift == 1) {
					bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
					bwBRCg.write(geneomeSeqID + "\t" + "vigor" + "\t");
					bw.write(frameshift == -1 ? "mRNA_with_minus_1_frameshift" : "mRNA_with_plus_1_frameshift" + "\t");
					bwBRCg.write(frameshift == -1 ? "mRNA_with_minus_1_frameshift" : "mRNA_with_plus_1_frameshift" + "\t");
					bw.write(Math.min(slippageBegin, slippageEnd) + "\t" + Math.max(slippageBegin, slippageEnd) + "\t");
					bwBRCg.write(Math.min(slippageBegin, slippageEnd) + "\t" + Math.max(slippageBegin, slippageEnd) + "\t");
					bw.write("." + "\t");
					bwBRCg.write("." + "\t");

					if (geneModel.getDirection().equals(Direction.FORWARD)) {
						bw.write("+" + "\t");
						bwBRCg.write("+" + "\t");
					} else {
						bw.write("-" + "\t");
						bwBRCg.write("-" + "\t");
					}
					bw.write("." + "\t");
					bwBRCg.write("." + "\t");
					idG = idGenerator.next();
					// bw.write(String.format("ID=%s;Parent=%s;",
					// idGenerator.next(), mRnaID));
					bw.write(String.format("ID=%s;Parent=%s;", idG, mRnaID));
					bwBRCg.write(String.format("ID=%s;Parent=%s;", idG, mRnaID));
					bw.write("\n");
					bwBRCg.write("\n");
				}
			}
		}
		bwBRCg.close();
	}
}
