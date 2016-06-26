package org.broadinstitute.cga.tools;

import java.io.Serializable;
import java.util.Comparator;
import java.util.Iterator;
import java.util.TreeSet;


public class MathSet extends TreeSet<Float> implements Serializable{
	public float sum = 0;
	public int n =0;
	public float max = Float.MIN_VALUE;
	public float min = Float.MAX_VALUE;
	
	public float getAverage() {return sum/((float) n);}
	public float getMean() {return getAverage();}
	
	static class FuzzyComparator implements Comparator<Float>, Serializable{

		/**
		 * 
		 */
		private static final long serialVersionUID = -7657114762143313468L;

		public int compare(Float o1, Float o2) {
			int comp = o1.compareTo(o2);
			if (comp == 0) return -1;
			return comp;
		}

	}

	public MathSet() {
		super(new FuzzyComparator());
	}



	@Override
	public synchronized boolean add(Float val) {
		float v = val;
		if (v > max) max = v;		
		if (v < min) min = v;
		n++;
		sum+=v;
		return super.add(val);
	}

	/*
	 * This returns the Av values in the form:
	 * .....
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString(){
		StringBuilder str = new StringBuilder();
		if (this.size() == 0) return "empty";
		float av = this.getAverage();
		str.append(String.valueOf(av));
		float stdDev = this.getStdDev();
		str.append('\t').append(stdDev);
		str.append('\t').append(stdDev / av);
		for (Float v: this) {
			str.append('\t');
			str.append(v.toString());
		}
		return str.toString();
	}

	
	/*
	 * Calculates Standard Deviation
	 */
	public float getStdDev(){
		float devSum = 0;
		float av = this.getAverage();
		float diff;
		for (Float d: this){
			diff = d - av;
			devSum += diff * diff;
		}
		return (float)Math.sqrt(devSum / (float)this.size());
	}
	
	
	
	private static final long serialVersionUID = 630911844067701546L;
	
	public float getMedian() {
		if (n == 0) throw new RuntimeException("getMedian called on empty MathVector");
		Iterator<Float> it = this.iterator();
		int half = n/2;
		float f  = Float.NaN;
		for ( int i =0; i <= half; i++ ){
			f = it.next();
		}
		
		return f;
	}

	public float getMax() {
		if (n == 0) throw new RuntimeException("getMax called on empty MathVector");
		return max;
	}
	public float getMin() {
		if (n == 0) throw new RuntimeException("getMin called on empty MathVector");
		return min;
	}
	
	public static void main (String [] args){
		try {
			

//			MathVector vm = new MathVector();
//			if (vm instanceof Serializable){
//				System.out.println("is serializable");
//			}else
//				System.out.println("is NOT serializable");
//
//			vm.add(9f);
//			vm.add(79f);
//			vm.add(20f);
//			vm.add(13f);
//			vm.add(224f);
//			vm.add(5f);
//			vm.add(66f);
//			vm.add(50f);
//			vm.add(4f);
//			vm.add(1f);
//			vm.add(11f);
//			System.out.println(vm.getMedian());
//			System.out.println(vm.size());
//			System.out.println(vm.getAverage());
		}catch (Exception e){
			e.printStackTrace();
		}
	}

	public float getLowerNthPercentile(float percentile) {
		if (n == 0) throw new RuntimeException("getNthPercentile called on empty MathVector");
		Iterator<Float> it = this.iterator();
		int cutoff = (int)( ((float)n) * percentile);
		float f  = Float.NaN;
		for ( int i =0; i <= cutoff; i++ ){
			f = it.next();
		}
		
		return f;
	}
	
	public Float getNthFromBottom(int n) {
		Iterator<Float> it = this.iterator();
		
		float f  = Float.NaN;
		for ( int i =0; i < n; i++ ){
			f = it.next();
		}
		
		return f;
	}
	
	public float getKurtosis(){
		float numSum = 0f;
		float denSum = 0f;
		float diff = 0;
		float mean = this.getMean();
		
		for (float f: this){
			diff = f - mean;
			denSum += (diff * diff );
			numSum += (diff * diff * diff * diff);			
		}
		
		float oneOverN = (1f / (float)this.size());
		return  (oneOverN * numSum) / (oneOverN*oneOverN*denSum*denSum) -3f;
		 
	
	}
}