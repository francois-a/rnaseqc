package org.broadinstitute.cga.rnaseq;

import org.broadinstitute.cga.tools.IOTools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by the Cancer Genome Analysis Group at the Broad Institute.
 * Author: David S. DeLuca
 * User: ddeluca
 * Date: 8/4/11
 * Time: 4:13 PM
 * Developed as part of the RNA-seq analysis efforts of the Broad Institute
 */
public class TranscriptListDoCResult {
    private Float averageCoverages;
    private Float averageCV;
    private Integer fiveEndCovCount;
    private Float avFiveEnd200;
    private Float avFiveEnd50;
    private Object avFiveEnd100;
    private Integer threeEndCovCount;
    private Float avThreeEnd200;
    private Float avThreeEnd50;
    private Float avThreeEnd100;
    private Integer numberOfGaps;
    private Integer cumulativeGapLength;
    private Float gapRatio;

    public Float getAverageCoverages() {
        return averageCoverages;
    }

    public Float getAverageCV() {
        return averageCV;
    }

    public Integer getFiveEndCovCount() {
        return fiveEndCovCount;
    }

    public Float getAvFiveEnd200() {
        return avFiveEnd200;
    }

    public Float getAvFiveEnd50() {
        return avFiveEnd50;
    }

    public Object getAvFiveEnd100() {
        return avFiveEnd100;
    }

    public Integer getThreeEndCovCount() {
        return threeEndCovCount;
    }

    public Float getAvThreeEnd200() {
        return avThreeEnd200;
    }

    public Float getAvThreeEnd50() {
        return avThreeEnd50;
    }

    public Float getAvThreeEnd100() {
        return avThreeEnd100;
    }

    public Integer getNumberOfGaps() {
        return numberOfGaps;
    }

    public Integer getCumulativeGapLength() {
        return cumulativeGapLength;
    }

    public Float getGapRatio() {
        return gapRatio;
    }

    public void loadSummaryFromFile(String doCTranscriptsResultsFile) throws IOException {
        ArrayList<String> values = IOTools.fileToList(doCTranscriptsResultsFile,null);
        this.averageCoverages = Float.valueOf(values.get(0));
        this.averageCV = Float.valueOf(values.get(1));
        this.fiveEndCovCount = Integer.valueOf(values.get(2));
        this.avFiveEnd200 = Float.valueOf(values.get(3));
        this.avFiveEnd50 = Float.valueOf(values.get(4));
        this.avFiveEnd100 = Float.valueOf(values.get(5));
        this.threeEndCovCount = Integer.valueOf(values.get(6));
        this.avThreeEnd200 = Float.valueOf(values.get(7));
        this.avThreeEnd50 = Float.valueOf(values.get(8));
        this.avThreeEnd100 = Float.valueOf(values.get(9));
        this.numberOfGaps = Integer.valueOf(values.get(10));
        this.cumulativeGapLength = Integer.valueOf(values.get(11));
        this.gapRatio = Float.valueOf(values.get(12));

    }
}
