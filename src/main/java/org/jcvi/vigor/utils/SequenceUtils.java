package org.jcvi.vigor.utils;
import org.jcvi.jillion.align.AminoAcidSubstitutionMatrix;
import org.jcvi.jillion.core.Sequence;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.nt.Nucleotide;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.Triplet;

import java.util.Iterator;
import java.util.regex.Pattern;
import java.util.stream.Stream;

public class SequenceUtils {


	public static Triplet getNextTriplet(Iterator<Nucleotide> iter){
		
		Nucleotide first= getNextNucleotide(iter);
		Nucleotide second= getNextNucleotide(iter);
		Nucleotide third= getNextNucleotide(iter);
		if(first==null || second ==null || third ==null){
			return null;
		}
		return Triplet.create(first, second, third);
	}
	
	private static Nucleotide getNextNucleotide(Iterator<Nucleotide> iter){
		if(!iter.hasNext()){
			return null;
		}
		Nucleotide n = iter.next();
		return n;
	}

	
	
	public static Iterator<Nucleotide> handleFrame(NucleotideSequence sequence, Frame frame) {
	    Iterator<Nucleotide> iter;
	    if(frame.onReverseStrand()){
	        iter = sequence.toBuilder().reverseComplement().iterator();
	          switch(frame){
                        case NEGATIVE_THREE:
                                        if(iter.hasNext()){
                                                iter.next();
                                        }
                        case NEGATIVE_TWO:
                                        if(iter.hasNext()){
                                                iter.next();
                                        }
                                        break;
                        default:
                                       
                                break;
                }
	    }else{
	        iter = sequence.iterator();
		switch(frame){
			case THREE:
					if(iter.hasNext()){
						iter.next();
					}
			case TWO:
					if(iter.hasNext()){
						iter.next();
					}
					break;
			default:
					break;
		}
	    }
	    return iter;
	}

	public static <T> Stream<String> steamOf(Sequence<T> sequence, int lineLength) {
		Pattern splitPattern = Pattern.compile(String.format("(?<=\\G.{%s})", lineLength));
		return splitPattern.splitAsStream(sequence.toString());
	}

	public static <T> String elipsedSequenceString(Sequence<T> sequence, int leading, int trailing) {
		String sequenceString = sequence.toString();
		if (leading + trailing >= sequence.getLength()) {
			return sequenceString;
		}
		return String.join("...",
				sequenceString.substring(0,leading),
				sequenceString.substring(sequenceString.length() -1 - trailing, sequenceString.length() - 1));
	}

	public static double computePercentSimilarity(ProteinSequence first, ProteinSequence second, int alignmentLength, AminoAcidSubstitutionMatrix matrix) {
		double similarity = 0;
		double maxLength = Math.max(first.getLength(), second.getLength());
		double minLength = Math.min(first.getLength(), second.getLength());
		// TODO gaps in the same place?
		double numberOfGaps = first.getNumberOfGaps() + second.getNumberOfGaps();
		double matches = 0;
		for (int i = 0; i < minLength; i++) {
			if (first.isGap(i) || second.isGap(i)) {
				continue;
			}
			matches += matrix.getValue(first.get(i), second.get(i)) > 0? 1: 0;
		}
		//similarity = (((int)(1000*(maxLength - numberOfGaps - misMatches))) / alignmentLength)/10;
        similarity = ((int)(1000*matches/ alignmentLength))/10;
		if(similarity>100) similarity=100;
        return similarity;
	}
}
