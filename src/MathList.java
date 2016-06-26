package org.broadinstitute.cga.tools;

import java.util.ArrayList;

public class MathList extends MathSet{
	/**
	 * 
	 */
	private static final long serialVersionUID = -4316111845894992137L;
	protected ArrayList<Float> vector = new ArrayList<Float>();
	
	@Override
	public synchronized boolean add(Float val) {
		vector.add(val);
		return super.add(val);
	}

	public Float get(int index){
		return vector.get(index);
	}

	
}

