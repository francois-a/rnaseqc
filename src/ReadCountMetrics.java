package org.broadinstitute.cga.rnaseq;

import org.broadinstitute.cga.picardbased.CountAligned;
import org.broadinstitute.cga.picardbased.CountAlignedExperimental;
import org.broadinstitute.cga.picardbased.LibraryComplexity;
import org.broadinstitute.cga.rnaseq.gatk.CountReadMetricsWalker;
import org.broadinstitute.cga.rnaseq.gatk.GATKTools;
import org.broadinstitute.cga.rnaseq.gatk.IntronicExpressionReadBlockWalker;
import org.broadinstitute.cga.tools.GCTFile;
import org.broadinstitute.cga.tools.IOTools;
import org.broadinstitute.cga.tools.Performance;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by the Cancer Genome Analysis Group at the Broad Institute.
 * Author: David S. DeLuca
 * User: ddeluca
 * Date: 8/2/11
 * Time: 11:23 AM
 * Developed as part of the RNA-seq analysis efforts of the Broad Institute
 */
public class ReadCountMetrics {

    private ArrayList<RNASeqMetrics.MetricSample> samples;
    private Object outDir;
    private String combinedGCTFileName;
    private float lowerExpressionCutoff;


    public ReadCountMetrics(ArrayList<RNASeqMetrics.MetricSample> samps, String outdir, float lowerExpressionCutoff) {
        this.samples = samps;
        this.outDir = outdir;
        this.combinedGCTFileName = outDir + "/genes.rpkm.gct";
        this.lowerExpressionCutoff = lowerExpressionCutoff;
    }


    /**
     * Performs the ReadWalker-based metrics. It also creates the RefGEn file for GATK.
     * Creates a rRNA interval file as well.
     * Runs readcounts on the rRNA intervals as well.
     *
     * Also runs a  routine to calculate library complexity and saves the results in
     * OUT_DIR+"/"+sampId+"/"+sampId+".libraryComplexity.txt"
     *
     * @param REF_GENOME
     * @param refGeneFile
     * @param  rRNAFile         If this file is a reference fasta, then bwa will run on it
     *                          If this file is a bam, then we will count the number of aligned reads
     *                          If this file is a *.list, we will run GATK and count reads in those intervals
     * @throws Exception
     */
    public void runReadCountMetrics(String REF_GENOME, String refGeneFile, String rRNAFile, boolean singleEnd, String bwa,
                                    int rRNAdSampleTarget, boolean skiprRNA, String gatkFlags, boolean strictMode) throws Exception {

        runRegionCounting(REF_GENOME, refGeneFile, gatkFlags, strictMode);

        if (!skiprRNA) {
            if (rRNAFile.endsWith(".list")) {
                countrRNAIntervals(REF_GENOME, rRNAFile, gatkFlags);
            } else {
                alignAndCountrRNA(rRNAFile, singleEnd, bwa, rRNAdSampleTarget);
            }
        } else {
            for (RNASeqMetrics.MetricSample samp: samples) {
                IOTools.stringToFile(samp.getrRNAMetricsFile(), "0\t0");
            }
        }
        calculateLibraryComplexity(singleEnd);
        combineMetricsPerSample();
        combineAcrossSamples();
    }

    /**
     * 1. combines the GCT Files into a single toplevel gct
     * 2. combines the metric results into a single file???? -> todo
     *
     * @return
     * @throws IOException
     */
    private String combineAcrossSamples() throws IOException {
        GCTFile merged = new GCTFile(samples.get(0).getGCTFile()); // get the first GCT File
        for (int i=1; i<samples.size(); i++) {
            GCTFile gct = new GCTFile(samples.get(i).getGCTFile());
            merged = merged.combinePreAligned(gct);
        }
        // fix sample ids
        for (int i=0; i<samples.size(); i++) {
            merged.setSampleId(i, samples.get(i).sampId);
        }
        // write result to file
        merged.toFile(this.getCombinedGCTFileName());
        return this.getCombinedGCTFileName();
    }


    private void combineMetricsPerSample() throws IOException {
        for (RNASeqMetrics.MetricSample samp: samples){
            combineMetrics(samp);
        }
    }


    /**
     * this concerns a metrics file for an individual sample. It combines the standard metrics outputted by the metrics walker
     * with rRNA metrics, library complexity
     *
     * @throws IOException
     */
    private void combineMetrics(RNASeqMetrics.MetricSample samp) throws IOException {

        // load RNA counts:
        String rRNACounts[] = IOTools.fileToString(samp.getrRNAMetricsFile()).split("\\t");
        int totalrRNA = Integer.valueOf(rRNACounts[0]);
        int uniquerRNA = Integer.valueOf(rRNACounts[1].trim());

        // load LC:
        String lc = LibraryComplexity.getFromFile(samp.getLibraryComplexityFile());

        // table 1
        BufferedReader in = new BufferedReader(new FileReader(samp.getTmpMetricsFile()));
        BufferedWriter out = new BufferedWriter(new FileWriter(samp.getMetricsFile()));
        out.write(in.readLine()); out.write("\tEstimated Library Size\n"); // header 1
        String data = in.readLine(); // capture the line containing total
        int totalReads = Integer.valueOf(data.split("\\t")[0]);
        out.write(data); out.write("\t"+lc +"\n"); // table 1 data, with injected library complexity

        // table 2:
        String header2 = in.readLine() + "\trRNA\trRNA rate"; //\tUnique non-rRNA\tUnique non-rRNA Rate"; // adding two fields to header
        out.write(header2); out.write('\n'); // header 2

        data = in.readLine();
        // int uniquMapped = Integer.parseInt(data.split("\\t")[2]);
        //int uniqueMappedNonrRNA = uniquMapped - uniquerRNA;

        data = data +"\t"+ totalrRNA+"\t"+((float)totalrRNA/(float)totalReads); //+"\t"+
        //uniqueMappedNonrRNA+"\t"+ ((float)uniqueMappedNonrRNA/(float)totalReads);
        out.write(data); out.write('\n');

        // copy the remainder of the tmp metrics file without modification
        String line = in.readLine();
        while (line!=null){
            out.write(line); out.write('\n');
            line = in.readLine();
        }
        out.close();
        in.close();
    }


    private void calculateLibraryComplexity(boolean singleEnd) throws IOException {
        for (RNASeqMetrics.MetricSample samp: samples) {
            System.out.println("Calculating library complexity for " + samp.sampId);
            Performance perf = new Performance("Libary Complexity Calculation Time");
            LibraryComplexity lc = new LibraryComplexity(!singleEnd);
            if (samp.hasList()) {
                ArrayList<String> files = IOTools.fileToList(samp.listFile, "\n");
                for (String f: files) {
                    lc.countFragmentsV2(f);
                }
            } else {
                lc.countFragmentsV2(samp.bamFile);
            }
            lc.toFile(samp.getLibraryComplexityFile());
            System.out.println(perf);
        }
    }

    /**
     *
     * @param rRNAFile
     */
    private void alignAndCountrRNA(String rRNAFile, boolean singleEnd, String bwa, int rRNAdSampleTarget) throws IOException, InterruptedException {
        for (RNASeqMetrics.MetricSample samp: samples) {

            Performance perf = new Performance("BWA based rRNA Estimation for " + samp.sampId, Performance.Resolution.minutes);
            System.out.println("Counting rRNA reads with BWA and " + rRNAFile);

            CountAlignedExperimental ca = new CountAlignedExperimental(rRNAFile,samp.getSampleDirectory());
            int totalReads = samp.getTotalReads();
            double rRNAdSampleRate = (double)rRNAdSampleTarget / (double)totalReads;

            ca.setdSampleRate(rRNAdSampleRate);
            if (bwa !=null) ca.BWA = bwa;
            if (samp.hasList()) {
                ArrayList<String> files = IOTools.fileToList(samp.listFile, "\n");
                for (String f: files){
                    ca.countBAM(f, singleEnd);
                }
            } else {
                ca.countBAM(samp.bamFile, singleEnd);
            }
            String result = ""+ ca.getTotalReads() + "\t"+ca.getUniqueReads();
            IOTools.stringToFile(samp.getrRNAMetricsFile(),result);
            System.out.println(perf);
        }
    }


    private void countrRNAIntervals(String refGenome, String rRNAIntervals, String gatkFlags) throws Exception {
        for (RNASeqMetrics.MetricSample samp: samples) {
            Performance perf = new Performance("CountReadMetricsWalker (rRNA) Runtime for "+samp.sampId, Performance.Resolution.minutes);
            System.out.println("Counting rRNA reads");
            // rRNAIntervals: file containing intervals
            GATKTools.runIntervalReadCounter(refGenome, samp.getBamFileOrList(), rRNAIntervals, samp.getrRNAMetricsFile(), gatkFlags);
            System.out.println(perf);
        }
    }


    private void runRegionCounting(String refGenome, String refGeneFile, String gatkFlags, boolean strictMode) throws Exception {
        for (RNASeqMetrics.MetricSample samp: samples){
            // replaced by "prepare for region counting ..
            Performance perf = new Performance("CountReadMetricsWalker Runtime", Performance.Resolution.minutes);
            String intervals = null;// "chr1:3530586-3534177";
            // previously: samp.getTmpMetricsFile(), which returned full path: <out_dir>/<sample_id>/<sample_id>.metrics.tmp.txt
            GATKTools.runIntronReadCount(refGenome, samp.getBamFileOrList(), null, refGeneFile, samp.getSampleDirectory()+"/"+samp.sampId, gatkFlags, strictMode);
            System.out.println(perf);
        }
    }


    public String getCombinedGCTFileName() {
        return combinedGCTFileName;
    }


    /**
     * Reads in the metrics result tables from the ReadCounter results that were prevsiouly run.
     * It reforms these tables into HTML and returns the resulting HTML
     *
     * @return
     * @throws IOException
     * @param metricTracker
     */
    public String getMetricsHTML(ArrayList<HashMap<String, String>> metricTracker) throws IOException {

        ArrayList<String[]> resultT1 = new ArrayList<String[]>();
        ArrayList<String[]> resultT2 = new ArrayList<String[]>();
        ArrayList<String[]> resultT3 = new ArrayList<String[]>();
        ArrayList<String[]> resultT4 = new ArrayList<String[]>();
        ArrayList<String[]> resultT5 = new ArrayList<String[]>();

        String[] header1 = null;
        String[] header2 = null;
        String[] header3 = null;
        String[] header4 = null;
        String[] header5 = null;

        StringBuilder metricTable = new StringBuilder();
        for (RNASeqMetrics.MetricSample s: samples){
            String metricFile = s.getMetricsFile();
            ArrayList<String[]> result= IOTools.fileToListofArrays(metricFile, "\\t");

            if (result == null){
                System.out.println("No Metrics data calculated. Skipping");
            } else {
                header1 = result.get(0);
                resultT1.add(result.get(1));
                header2 = result.get(2) ;
                resultT2.add(result.get(3));
                header3 = result.get(4);
                resultT3.add(result.get(5));
                header4 = result.get(6);
                resultT4.add(result.get(7));
                header5 = result.get(8);
                resultT5.add(result.get(9));
            }
        }

        if (header1 != null){
            metricTable.append("\n<h3>Total Reads</h3>\n");
            metricTable.append(RNASeqMetrics.getMetricTable(samples, resultT1, header1, CountReadMetricsWalker.TABLE1_FORMAT + "0", metricTracker)); // total reads and duplication
            metricTable.append("\n<p><b>Total Purity Filtered Reads Sequenced</b> are filtered for vendor fail flags and exclude alternative alignment reads. <b>Alternative Aligments</b> are duplicate read entries " +
                    "providing alternative coordinates. <b>Failed Vendor QC Check</b> are reads which have been designated as failed by the sequencer. <b>Read Length</b>" +
                    " is the maximum length found for all reads. <b>Estimated Library Size</b> is the number of " +
                    "expected fragments based upon the total number of reads and duplication rate assuming a Poisson distribution.</p>\n");
            metricTable.append("\n<h3>Mapped Reads</h3>\n");
            metricTable.append(RNASeqMetrics.getMetricTable(samples, resultT2, header2, CountReadMetricsWalker.TABLE2_FORMAT+"01", metricTracker));
            metricTable.append("\n<p><b>Mapped</b> reads are those that were aligned. <b>Mapping Rate</b> is per total reads." +
                    " <b>Mapped Unique</b> are both aligned as well as non-duplicate reads. <b>Mapped Unique Rate of Total</b> is" +
                    " per total reads. <b>Unique Rate of Mapped</b> are unique reads divided by all mapped reads. <b>Duplication Rate of Mapped</b> is the duplicate read divided by total mapped reads. " +
                    " <b>Base Mismatch Rate</b> is the number of bases not matching the reference divided by the total number of aligned bases. " +
                    "<b>rRNA</b> reads are non-duplicate and duplicate reads aligning to rRNA regions " +
                    "as defined in the transcript model definition. <b>rRNA Rate</b> is per total reads. </p>\n");

            metricTable.append("\n<h3>Mate Pairs</h3>\n");
            metricTable.append(RNASeqMetrics.getMetricTable(samples, resultT3, header3, CountReadMetricsWalker.TABLE3_FORMAT , metricTracker));
            metricTable.append("\n<p>" +
                    //"Mapped Pairs\tEnd 1 Mapping Rate\tEnd 2 Mapping Rate\tEnd 1 Mismatch Rate\tEnd 2 Mismatch Rate\tFragment Length Mean\tFragment Length StdDev\tChimeric Pairs");
                    "<b>Mapped Pairs</b> is the total number of pairs for which both ends map. <b>Unpaired Reads</b> are the number of reads that are lacking a mate. <b>End 1/2 Mapping Rate</b> is the number of mapped divided by the total number of End1/End2 reads." +
                    " <b>End 1/2 Mismatch Rate</b> is the number of End 1 and 2 bases not matching the reference divided by the total number of mapped End 1 and 2 bases. "+
                    "<b>Fragment Length Mean/StdDev</b> is the mean distance, standard deviation between the start of an upstream read and the end of the downstream one. Only fragments contained within single exons are used. "+
                    "<b>Chimeric Pairs</b> are pairs whose mates map to different genes. " +
                    "</p>\n");

            metricTable.append("\n<h3>Transcript-associated Reads</h3>\n");
            metricTable.append(RNASeqMetrics.getMetricTable(samples, resultT4, header4, CountReadMetricsWalker.TABLE4_FORMAT , metricTracker));
            metricTable.append("\n<p>All of the above rates are per mapped read. <b>Intragenic Rate</b> refers to the fraction of reads " +
                    " that map within genes (within introns or exons). <b>Exonic Rate</b> is the fraction mapping within exons. <b>Intronic Rate</b> is the fraction " +
                    "mapping within introns. <b>Intergenic Rate</b> is the fraction mapping in the genomic space between genes." +
                    " <b>Split Reads</b> is the number of reads spanning an exon exon junction. " +
                    " <b>Expression Profile Efficiency</b> is the ratio of exon reads to total reads. <b>Transcripts/Genes Detected </b> is the number of " +
                    "transcripts/Genes with at least "+ IntronicExpressionReadBlockWalker.LOWER_READ_COUNT_CUTOFF+" reads. </p>\n");

            metricTable.append("\n<h3>Strand Specificity</h3>\n");
            metricTable.append( RNASeqMetrics.getMetricTable(samples, resultT5, header5, CountReadMetricsWalker.TABLE5_FORMAT, metricTracker));
            metricTable.append("\n<p><b>End 1/2 Sense</b> are the number of End 1 or 2 reads that were sequenced in the sense direction. " +
                    "Similarly, <b>End 1/2 Antisense</b> are the number of End 1 or 2 reads that were sequenced in the antisense direction."+
                    "<b>End 1/2 Sense %</b> are percentages of intragenic End 1/2 reads that were sequenced in the sense direction. </p>\n");

            IOTools.stringToFile(this.outDir+"/countMetrics.html", metricTable.toString()); // just ouput it so we can peek
        } else {
            metricTable.append("none calculated");
        }
        return metricTable.toString();
    }

}
