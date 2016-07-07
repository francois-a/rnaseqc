package org.broadinstitute.cga.rnaseq;

import org.broadinstitute.cga.tools.ObjectCounter;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: ddeluca
 * Date: 3/25/11
 * Time: 1:57 PM
 */
public class DoCMetrics {

    private int totalCoverage;
    private int length;
    private float meanCoverage;

    private boolean has3EndCov;
    private boolean has5EndCov;

    private int cumulativeGapLength;
    private int numberOfGaps;
    private int threeEnd100;
    private int threeEnd50;
    private int threeEnd10;

    private int fiveEnd100;
    private int fiveEnd50;
    private int fiveEnd10;
    private float stdDev;
    private ArrayList<Integer> coverage;
    
    private ObjectCounter<Integer> gapHistogram = null;
    private static  int gapHistogramGreaterThanThresh = 100;
    private int gapHistogramGreaterThanCount = 0;

    public void setFiveEnd400(float fiveEnd400) {
        this.fiveEnd400 = fiveEnd400;
    }

    public void setThreeEnd400(float threeEnd400) {
        this.threeEnd400 = threeEnd400;
    }

    private float fiveEnd400;
    private float threeEnd400;


    public DoCMetrics(TranscriptDoCResult result) {
        set5End10Cov(result.get5EndCoverage(10));
        set5End50Cov(result.get5EndCoverage(50));
        set5End100Cov(result.get5EndCoverage(100));
        setFiveEnd400(result.get5EndCoverage(400));

        set3End10Cov(result.get3EndCoverage(10));
        set3End50Cov(result.get3EndCoverage(50));
        set3End100Cov(result.get3EndCoverage(100));
        setThreeEnd400(result.get3EndCoverage(400));

        setNumberOfGaps(result.getNumberOfGaps());
        setCumulativeGapLength(result.getCumlativeGapLength());

        setHas5EndCov(result.has5EndCoverage(100));
        setHas3EndCov(result.has3EndCoverage(100));

        setTotalCoverage(result.getTotalCoverage());
        setLength(result.getSize());
        setStdDev(result.getStdDev());

        meanCoverage = (float)this.totalCoverage / (float) this.length;

        this.coverage = result.get5PrimeTo3PrimeCoverage();
        
        createGapHistogram(result);
    }

    private void createGapHistogram(TranscriptDoCResult result) {
        this.gapHistogram = new ObjectCounter<Integer>();
    }

    public boolean has3EndCov() {
        return has3EndCov;
    }

    public boolean has5EndCov() {
        return has5EndCov;
    }

    public void setCumulativeGapLength(int cumulativeGapLength) {
        this.cumulativeGapLength = cumulativeGapLength;
    }

    public void setNumberOfGaps(int numberOfGaps) {
        this.numberOfGaps = numberOfGaps;
    }

    public void set3End100Cov(int cov) {
        this.threeEnd100 = cov;
    }

    public void set3End50Cov(int cov) {
        this.threeEnd50 = cov;
    }

    public void set3End10Cov(int cov) {
        this.threeEnd10 = cov;
    }

    public void set5End100Cov(int cov) {
        this.fiveEnd100 = cov;
    }

    public void set5End50Cov(int cov) {
        this.fiveEnd50 = cov;
    }

    public void set5End10Cov(int cov) {
        this.fiveEnd10 = cov;
    }

    public int getNumberOfGaps() {
        return this.numberOfGaps;
    }
    
    public int getCumulativeGapLength() {
        return this.cumulativeGapLength;
    }

    /**
     * This is the ratio between the cumulative gap lengthg and transcript length
     * @return
     */
    public float getGapRatio() {
        return (float)this.cumulativeGapLength / (float)this.length;
    }

    public float getThreeEnd100Norm() {
        return this.threeEnd100 / this.getMeanCoverage();
    }

    public float getThreeEnd50Norm() {
        return this.threeEnd50 / this.getMeanCoverage();
    }

    public float getThreeEnd10Norm() {
        return this.threeEnd10 / this.getMeanCoverage();
    }

    public float getFiveEnd100Norm() {
        return this.fiveEnd100 / this.getMeanCoverage();
    }

    public float getFiveEnd50Norm() {
        return this.fiveEnd50 / this.getMeanCoverage();
    }

    public float getFiveEnd10Norm() {
        return this.fiveEnd10 / this.getMeanCoverage();
    }

    public float getMeanCoverage() {
        return meanCoverage;
    }

    public void setHas3EndCov(boolean endCoverage) {
        this.has3EndCov = endCoverage;
    }

    public void setHas5EndCov(boolean endCoverage) {
        this.has5EndCov = endCoverage;

    }

    public float getCV(){
        float std = this.getStdDev();
        float mean = this.getMeanCoverage();
        return std/mean;
    }

    public float getStdDev(){
        return this.stdDev;
    }

    public int getTotalCoverage() {
        return totalCoverage;
    }

    public void setTotalCoverage(int totalCoverage) {
        this.totalCoverage = totalCoverage;
    }

    public void setLength(int length) {
        this.length = length;
    }

    public void setStdDev(float stdDev) {
        this.stdDev = stdDev;
    }

    public int getLength() {
        return length;
    }

    public boolean isHas3EndCov() {
        return has3EndCov;
    }

    public boolean isHas5EndCov() {
        return has5EndCov;
    }

    public int getThreeEnd100() {
        return threeEnd100;
    }

    public int getThreeEnd50() {
        return threeEnd50;
    }

    public int getThreeEnd10() {
        return threeEnd10;
    }

    public int getFiveEnd100() {
        return fiveEnd100;
    }

    public int getFiveEnd50() {
        return fiveEnd50;
    }

    public int getFiveEnd10() {
        return fiveEnd10;
    }

    public ArrayList<Integer> getCoverage() {
        return coverage;
    }

    public float getFiveEnd400() {
        return fiveEnd400;
    }

    public float getThreeEnd400() {
        return threeEnd400;
    }
}
