package org.jcvi.vigor.component;


import java.util.Comparator;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.internal.core.util.JillionUtil;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data

public class Exon implements Cloneable {

	private Range range;
	private Frame frame;
	private AlignmentFragment alignmentFragment;
	private boolean is_5p_adjusted=false;
	private boolean is_3p_adjusted=false;
	private String insertionString;
	private String replacementString;
	private Frame sequenceFrame;
	
	public Exon(){
		
	}
	public Exon(Range range,Frame frame){
		this.range = range;
		this.frame = frame;
	}
	public Exon clone(){
		
		Exon exon=null;
		try {
			exon = (Exon)(super.clone());
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return exon;
	}
	/*
	 public static Exon deepClone(Exon exon) {
		   try {
		     ByteArrayOutputStream baos = new ByteArrayOutputStream();
		     ObjectOutputStream oos = new ObjectOutputStream(baos);
		     oos.writeObject(exon);
		     ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
		     ObjectInputStream ois = new ObjectInputStream(bais);
		     return (Exon)(ois.readObject());
		   }
		   catch (Exception e) {
		     e.printStackTrace();
		     return null;
		   }
	 }*/
	 
	public enum Comparators implements Comparator<Exon>{
		Descending{
			@Override
			  public int compare(Exon e1, Exon e2) {
                return -1 * JillionUtil.compare(e1.getRange().getBegin(), e2.getRange().getBegin());
            }
		}
	,
	
	Ascending{

        @Override
        public int compare(Exon e1, Exon e2) {
            return JillionUtil.compare(e1.getRange().getBegin(),e2.getRange().getBegin());
        }
        
    };
	 
	}

		
}
