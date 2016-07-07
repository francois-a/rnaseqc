package org.broadinstitute.cga.rnaseq;

import org.broadinstitute.cga.tools.MathSet;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: ddeluca
 * Date: 3/22/11
 * Time: 1:25 PM
 *
 * This class helps organize the results of the DoC analysis as output by GATK.
 * Each line in the result file corresponds to a base and it's DoC data.
 *
 * It handles some of the metrics associated with the per-base DoC.:
 *          End coverage
 *          Gaps
 *
 */
public class TranscriptDoCResult {
    static int MIN_GAP_SIZE = 5;
    static int MAX_COV_FOR_GAP = 0;

    String id, strand;
    ArrayList<Integer> values = new ArrayList<Integer>();
    ArrayList<Integer> gapLengths = null;

    TranscriptDoCResult(String id, String strand) {
        this.id = id;
        this.strand = strand;
    }

    void writeToFile(String outDir) throws IOException {
        BufferedWriter out = new BufferedWriter(new FileWriter(outDir+"/"+this.id+".DoC"));
        for (Integer val: values) {
            out.write(val.toString());
            out.write('\n');
        }
        out.close();
    }

    /**
     * line has the form:
     * transcript id<tab>gene name<tab> Locus   Total_Depth     Average_Depth_sample    Depth_for_an5_27th
     * @param line
     */
    void addLine(String line) {
        values.add(new Integer(line.split("\\t")[3]));
    }

    /**
     * calculates the average per base coverage for 'length' bases at the end
     *
     * @param length
     * @return
     */

    public int get5EndCoverage(int length) {
        if (this.strand.equals("1")|| this.strand.equals("+")) {
            // positive strand, so our 5 prime end is the left end
            return getLeftEndCoverage(length);
        } else {
            // negative strand, so our 5 prime end is right end
            return getRightEndCoverage(length);
        }
    }


    public ArrayList<Integer> get5PrimeTo3PrimeCoverage(){
        if (this.strand.equals("1")|| this.strand.equals("+")) {
            // positive strand, so our 5 prime end is the left end
            return values;
        } else {
            // negative strand, so our 5 prime end is right end
            int l = values.size();
            ArrayList<Integer> cov = new ArrayList<Integer>(l);
            for (int i = 0 ; i < l; i++) {
                int backDex = l-1-i;
                cov.add(values.get(backDex));
            }
            return cov;
        }
    }
    /**
     * calculates the average per base coverage for 'length' bases at the end
     *
     * @param length
     * @return
     */
    public int get3EndCoverage(int length) {
        if (this.strand.equals("1") || this.strand.equals("+")) {
            // positive strand, so our 5 prime end is the left end
            return getRightEndCoverage(length);
        } else {
            // negative strand, so our 5 prime end is right end
            return getLeftEndCoverage(length);
        }
    }

    int getLeftEndCoverage(int length){
        int sum = 0;
        if (values.size() < length) {
            return -1;
        }
        for (int i=0 ; i < length; i++) {
            sum+= values.get(i);
        }
        return sum / length;
    }

    int getRightEndCoverage(int length) {
        int sum = 0;
        if (values.size() < length) {
            return -1;
        }
        for (int i=1 ; i <= length; i++) {
            sum+= values.get(values.size() - i);
        }
        return sum / length;
    }

    boolean has5EndCoverage(int length) {
        return get5EndCoverage(length) > 0;
    }

    boolean has3EndCoverage(int length) {
        return get3EndCoverage(length) > 0;
    }

    void calcGapLengths() {
        this.gapLengths = new ArrayList<Integer>();
        int thisGap =-1;
        for (int i = 0; i < values.size(); i++){
            if (values.get(i) <= MAX_COV_FOR_GAP){
                // we are in a gap
                if (thisGap == -1){
                    // we are starting a new gap
                    thisGap =0;
                }
                thisGap++; // extend the gap
            } else {
                // no gap.
                // if we are ending a gap:
                if (thisGap > 0) {
                    // save last gap
                    if (thisGap >= MIN_GAP_SIZE) {
                        this.gapLengths.add(thisGap);
                    }    
                    thisGap = -1;
                }
            }
        }
        // check whether the end is a gap
        if (thisGap > 0) {
            if (thisGap >= MIN_GAP_SIZE) this.gapLengths.add(thisGap);
        }
    }

    public int getNumberOfGaps() {
        if (this.gapLengths == null) calcGapLengths();
        return this.gapLengths.size();
    }

    /**
     *
     * @return the cumulative length of the gaps in this tanscript
     */
    public int getCumlativeGapLength() {
        if (this.gapLengths == null) calcGapLengths();
        int sum = 0;
        for (Integer l: this.gapLengths){
            sum = sum + l;
        }
        return sum;
    }

    public int getTotalCoverage() {
        int tot = 0;
        for (Integer val: values) {
            tot+=val;
        }
        return tot;
    }

    public int getSize() {
        return values.size();
    }

    public float getStdDev() {
        MathSet m = new MathSet();
        for (float f: values) {
            m.add(new Float(f));
        }
        return m.getStdDev();
    }

}

