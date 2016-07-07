package org.broadinstitute.cga.rnaseq;

import java.awt.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.*;

import org.broadinstitute.cga.rnaseq.gatk.GATKTools;
import org.broadinstitute.cga.tools.GCTFile;
import org.broadinstitute.cga.tools.IOTools;
import org.broadinstitute.cga.tools.ObjectCounter;
import org.broadinstitute.cga.tools.Performance;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;


/**
 *  This class encapulates the job of running GATK depth of coverage analysis on a per-base basis.
 */
public class PerBaseDoC {

	static String ENSEMBL_BASEURL = "http://uswest.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=";
    static int MAX_GNUPLOTS = 100;
	static int DISTANCE_FROM_3PRIME_PLOT_BUCKETS = 2000;
    static int GAP_LENGTH_HISTOGRAM_BUCKETS = 500;
	static String USAGE = "java -jar GtexDoc.jar [--noGATK] No.OfGenes EXPRESSION_FILE GENOME_MAP REF_GENOME BAM_FILE OUT_DIR";
    private ArrayList<RNASeqMetrics.MetricSample> samples;
    static  String EXPR_SUFFIX = ".transcripts.list";
	String sampleId;
    boolean noGatk = false;
    private String outDir;
    enum ExpressionLevel {low, medium, high};  //, all};  // breaks pre-1.1.9 version


    public PerBaseDoC(ArrayList<RNASeqMetrics.MetricSample> bams, String OUT_DIR) {
        this.outDir = OUT_DIR;
        this.samples = bams;
    }


    public static void main(String[] args) {
        try {
            if (args.length <7) {
                System.out.println("USAGE:\t"+USAGE);
                return;
            }

            // GET ARGUMENTS****************************
            int a = 0;
            boolean noGatk = false;
            if (args[0].equals("--noGATK")){
                noGatk = true;
                a++;
            }
            int MAX = Integer.parseInt(args[a++]);
            String EXPRESSION_FILE = args[a++];
            String GENOME_MAP = args[a++];
            String REF_GENOME = args[a++];
            String BAM_FILE = args[a++];
            String OUT_DIR = args[a++];
            // ******************************************

            String sampleId = EXPRESSION_FILE.replace(".geneexpression.txt", "");
            if (sampleId.equals(EXPRESSION_FILE)) {
                sampleId="AnonymousSample";
            }

/*
            ArrayList<String> bams  = new ArrayList<String>();
            bams.add(BAM_FILE);
            PerBaseDoC doc = new PerBaseDoC(bams, OUT_DIR);
            doc.createExpressionStratifiedTranscriptLists(MAX, theseSamplesGCT,LOWER_EXPR_CUTOFF);
            doc.prepareIntervals(TRANSCRIPT_MODEL,transcriptTypeField);

            doc.runDoC(REF_GENOME, dSampleTarget, gatkFlags);
            // Create Index.html

            doc.createReports(TRANSCRIPT_MODEL,transcriptTypeField, details, this.gcFile, endLength);


            createTopLevelIndex(bams,MAX, metricTable, corrTable,endLength, OUT_DIR+"/index.html");
            System.out.println("Done.");
            */
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("USAGE:\t"+USAGE);
        }
    }


    /**
     * Generates the data for the visual coverage by distance from 3' end of transcript.
     * Saves the data in a file : OUT_DIR/meanCovByPosition.txt
     *      ( since this is PerBAseDoC, OUT_DIR will typically be the directory for the specific sample)
     *
     * The resulting file contains one value per line. The values are the average coverage at that position. The position
     * order is from 3' to 5' prime
     * @throws IOException
     */
    private void calculateAverageCoverageByPosition(TranscriptList transcripts, String outfile) throws IOException{
        // create a 1000 buckets
        int buckets = DISTANCE_FROM_3PRIME_PLOT_BUCKETS;

        int[] counts = new int[buckets];
        int[] covs = new int[buckets];

        Arrays.fill(counts, 0);
        Arrays.fill(covs,0);

        for (Transcript t: transcripts) {
            ArrayList<Integer> cov = t.getCoverage();
            for (int i = 0 ; i < buckets&& i< cov.size(); i++) {
                int backDex = cov.size()-i-1;
                counts[i]++;
                covs[i]+=cov.get(backDex);  // the resulting values are from 3' end to the 1000 positions from 3'
            }
        }

        float[] averageCovs = new float[buckets];
        for (int i = 0 ; i < buckets; i++) {
            averageCovs[i] = (float)covs[i] / (float)counts[i];
        }

        //save to file:
        BufferedWriter out = new BufferedWriter (new FileWriter(outfile));
        for (int i = 0 ; i < buckets ; i++) {
            out.write("" + averageCovs[i] + "\n");
        }
        out.close();
    }


    /**
     * Generates the data for the visual coverage by distance from 3' end of transcript.
     * Saves the data in a file : OUT_DIR/meanCovByPosition.txt
     *      ( since this is PerBAseDoC, OUT_DIR will typically be the directory for the specific sample)
     *
     * The resulting file contains one value per line. The values are the average coverage at that position. The position
     * order is from 3' to 5' prime
     * @throws IOException
     */
    private void calculateAverageCoverageByNormalizedPosition(TranscriptList transcripts, String outfile) throws IOException{
        int buckets = 100;

        int[] counts = new int[buckets];
        int[] covs = new int[buckets];

        Arrays.fill(counts, 0);
        Arrays.fill(covs,0);

        for (Transcript t: transcripts) {
            if (t.getCoverage().size()>0){
                ArrayList<Integer>cov = scale(t.getCoverage(), buckets);
                for (int i = 0 ; i < cov.size(); i++) {
                    counts[i]++;
                    covs[i]+=cov.get(i);
                }
            }
        }

        float[] averageCovs = new float[buckets];
        for (int i = 0 ; i < buckets; i++) {
            if(counts[i] == 0) {
                averageCovs[i] = 0;
            } else {
                averageCovs[i] = (float)covs[i] / (float)counts[i];
            }
        }

        //save to file:
        BufferedWriter out = new BufferedWriter (new FileWriter(outfile));
        for (int i = 0 ; i < buckets ; i++) {
            out.write("" + averageCovs[i] + "\n");
        }
        out.close();
    }

    
    private static ArrayList<Integer> scale(ArrayList<Integer> original, int newLength) {
        ArrayList<Integer> scaled = new ArrayList<Integer>(newLength);

        float scaleFactor = (float)original.size() / (float)newLength;
        for (int i = 0 ; i < newLength; i++) {
            int start = (int)(i*scaleFactor);
            int stop = (int)((i+1)*scaleFactor);

            if (start == stop) {
                stop = start+1;
            }
            
            int sum = 0;
            int n = 0;
            for (int j = start; (j< stop) && (j < original.size()); j++) {
                sum += original.get(j);
                n++;
            }
            if (n ==0) throw new RuntimeException("Could not scale vector from length " + original.size() + " to " + newLength + " with scale factor: " + scaleFactor);

            scaled.add(sum/n);
        }

        return scaled;
    }


    /**
     * Generates the data for the visual coverage by distance from 3' end of transcript.
     * Saves the data in a file : OUT_DIR/meanCovByPosition.txt
     *      ( since this is PerBAseDoC, OUT_DIR will typically be the directory for the specific sample)
     *
     * The resulting file contains one value per line. The values are the average coverage at that position. The position
     * order is from 3' to 5' prime
     * @throws IOException
     */
    private void calculateGapLengthHistogram(TranscriptList transcripts, String outfile) throws IOException{

        ObjectCounter<Integer> gapHist = new ObjectCounter<Integer>();
        int topBin = 0;
        int topBinLength = GAP_LENGTH_HISTOGRAM_BUCKETS;

        for (Transcript t: transcripts) {
            topBin += t.binGaps(gapHist,topBinLength);
        }

        //save to file:
        BufferedWriter out = new BufferedWriter (new FileWriter(outfile));
        for (int i = 1 ; i < topBinLength ; i++) {
            out.write("" + gapHist.getCount(new Integer(i)) + "\n");
        }
        out.write("" + topBin + "\n");
        out.close();
    }


    public static void cleanDir(String oUT_DIR) {
        File dir = new File(oUT_DIR);
        dir.mkdir();
        for (String file: dir.list()){
            if (file.endsWith(".DoC")){
                file = dir.getAbsolutePath() + "/" + file;
            }
        }
        System.out.println("Writing DoC per gene into: "+dir.getAbsolutePath());
    }


    /**
     * Creates one big long HTML page for all the plots, as well as an individual page for each transcript
     *
     * @param outDir
     * @param genes
     * @throws IOException
     */
    private void createDetailedHTML(String outDir, Collection<Transcript> genes, boolean useGC) throws IOException{
        BufferedWriter outPlots = new BufferedWriter(new FileWriter(outDir+"/plots.html"));

        outPlots.write("<html>\n<head><title>Depth of Coverage Report for "+this.sampleId+"</title>"+RNASeqMetrics.getStyle()+"</head>\n");
        outPlots.write("<body>\n");

        for (Transcript g: genes) {
            String entryHTML = createDetailedHTML(g, useGC);
            outPlots.write(entryHTML);

            BufferedWriter out = new BufferedWriter(new FileWriter(outDir+"/"+g.getTranscriptId()+".html"));
            out.write("<html>\n<head><title>Depth of Coverage Report for "+g.getTranscriptId()+"</title>"+RNASeqMetrics.getStyle()+"</head>\n");
            out.write("<body>\n");
            out.write(entryHTML);
            out.write("<h5>"+new Date()+"</h5>\n");
            out.write("</body>\n</html>");
            out.close();
        }

        outPlots.write("<h5>"+new Date()+"</h5>\n");
        outPlots.write("</body\n</html>");
        outPlots.close();
    }

    /**
     * Formats an report entry in HTML. This corresponds to a single transcript
     * Does not include body or html header
     *
     * @param g
     * @return
     * @throws IOException
     */
    private static String createDetailedHTML(Transcript g, boolean useGC) throws IOException{

        DoCMetrics m = g.getDoCMetrics();
        StringBuilder str = new StringBuilder();

        // ***************************************************************
        // TRANSCRIPT AND COVERAGE (expression file stuff)
        str.append("<h2>").append(g.getAnnotation().getName()).append(" / ").append(g.getTranscriptId()).append("</h2></a>\n");
        str.append("<h3>Expression file info</h3>\n");
        str.append("<table>");
        str.append("<tr><td>Transcript Id:</td><td><a target='_new' href='"+ENSEMBL_BASEURL + g.getTranscriptId() + "'>" + g.getTranscriptId() + "</a></td></tr>");
        str.append("<tr><td>Transcript Length:</td><td align='right'>");

        str.append(g.getLength());
        str.append("</td></tr><tr><td>Total Coverage:</td><td align='right'>");
        str.append(m.getTotalCoverage());
        str.append("</td></tr><tr><td>Average Coverage:</td><td align='right'>");
        str.append(m.getMeanCoverage());
        str.append("</td></tr><tr><td>Standard Deviation:</td><td align='right'>");
        str.append("" + String.format("%.2f", m.getStdDev()));
        str.append("</td></tr><tr><td>Coefficient of Variation:</td><td align='right'>");
        str.append(""+m.getCV());

        str.append("</td>\n");

        str.append("</tr></table>\n");

        // ***************************************************************
        // ENDS
        str.append("<h3>Ends</h3>\n");
        str.append("<table border =1> <tr> <th>End</th><th>10 Bases</th><th>50 Bases</th><th>100 Bases</th>");
        str.append("<th>10 Base Norm</th><th>50 Base Norm</th><th>100 Base Norm</th></tr>\n");
        str.append("<tr><td>5'</td>");
        str.append("<td align='right'>"+m.getFiveEnd10() + "</td>");
        str.append("<td align='right'>"+m.getFiveEnd50() + "</td>");
        str.append("<td align='right'>").append(m.getFiveEnd100()).append("</td>");
        str.append("<td align='right'>"+ String.format("%.3f",(m.getFiveEnd10()/ (m.getMeanCoverage()))) + "</td>");
        str.append("<td align='right'>"+String.format("%.3f",(m.getFiveEnd50()/(m.getMeanCoverage()))) + "</td>");
        str.append("<td align='right'>"+String.format("%.3f",(m.getFiveEnd100()/(m.getMeanCoverage()))) + "</td>");

        str.append("</tr>\n<tr>\n<td>3'</td>");

        str.append("<td align='right'>"+m.getThreeEnd10() + "</td>");
        str.append("<td align='right'>").append(m.getThreeEnd50()).append("</td>");
        str.append("<td align='right'>"+m.getThreeEnd100() + "</td>");
        str.append("<td align='right'>"+ String.format("%.3f",(m.getThreeEnd10()/ (m.getMeanCoverage()))) + "</td>");
        str.append("<td align='right'>"+String.format("%.3f",(m.getThreeEnd50()/(m.getMeanCoverage()))) + "</td>");
        str.append("<td align='right'>"+String.format("%.3f",(m.getThreeEnd100()/(m.getMeanCoverage()))) + "</td>");
        str.append("</tr>\n</table>\n");

        str.append("<p>End coverage is averaged per base. <i>Base Norm</i> is normalized by dividing by the transcript's average coverage.\n</p>\n");

        // ***************************************************************
        // GAPS
        str.append("<h3>Gaps</h3>\n");

        str.append("<table>");
        str.append("<tr><td>Number of Gaps:</td><td align='right'>"+m.getNumberOfGaps()+"</td></tr>");
        str.append("<tr><td>Cumulative Gap Length:</td><td align='right'>");
        str.append(""+m.getCumulativeGapLength());
        str.append("</td></tr><tr><td>Gap Percentage:</td><td align='right'>");
        str.append(""+String.format("%.1f",(m.getGapRatio() * 100f))+"%");
        str.append("</td></tr><tr><td>Average Coverage:</td><td align='right'>");
        str.append(String.format("%.2f",m.getMeanCoverage()));
        str.append("</td>\n");
        str.append("</tr></table>\n");

        str.append("<p>A gap is defined as a stretch of at least " +TranscriptDoCResult.MIN_GAP_SIZE+" bases having a coverage of no more than " + TranscriptDoCResult.MAX_COV_FOR_GAP+ "\n</p>\n");

        // ***************************************************************
        // OTHER STUFF (intervals, geneId, biotype, gc content, sequence	
        str.append("<h3>Genomic Intervals</h3>"+g.getExonList() + "<br>\n");

        str.append("<h3>Other Annotations</h3>\n");
        str.append("<table><tr><td>Gene ID:</td><td>");
        str.append("<a target='_new' href='"+ENSEMBL_BASEURL + g.getAnnotation().getGeneId() + "'>" + g.getAnnotation().getGeneId() + "</a>");
        str.append("</td></tr><tr><td>Strand:</td><td>");
        str.append(g.getAnnotation().getStrand());
        str.append("</td></tr><tr><td>Type:</td><td>");
        str.append(g.getAnnotation().getType());
        str.append("</td></tr>");
        if(useGC) {
            str.append("<tr><td>GC Content:</td><td>");
            str.append(g.getAnnotation().getGCContent());
            str.append("</td></tr>");
        }
        str.append("</table>\n");
        if((new File(g.getTranscriptId()+".png")).exists()) {
            str.append("<img src = '"+g.getTranscriptId()+".png'><br>");
            str.append("<br>\n");
        }
        // sequence
        str.append("<p>\n");
        str.append(g.getAnnotation().getSequence());
        str.append("\n</p>\n");

        return str.toString();
    }



    /**
     * Creates a png plot using GNUPlot
     * It opens gnuplot as a process and steams the commands to it
     *
     * Because of the storage requirements for creating so many plots, the number is restricted to MAX_GNUPLOTS
     * @throws IOException
     * @throws InterruptedException
     */
    private void createGNUPlots(TranscriptList transcripts, String sampleDir) throws IOException, InterruptedException {
        File dir = new File(sampleDir);
        ProcessBuilder pb = new ProcessBuilder("gnuplot","-persist");
        pb.redirectErrorStream(true);
        Process pro = pb.start();

        // create thread to read out the standard out
        IOTools.StdOGobbler gobbler = new IOTools.StdOGobbler(pro);
        gobbler.start();

        // pipe commands into process:
        OutputStreamWriter out = new OutputStreamWriter(pro.getOutputStream());
        ArrayList<String> tempFiles = new ArrayList<String> (MAX_GNUPLOTS);
        int count = 0;
        for (Transcript gene: transcripts){
            String f = gene.getTranscriptId() + ".DoC";
            // make it output to a png:
            out.write("set terminal png;\n");

            // set axis labels
            out.write("set xlabel \"Position in Gene\"\n");
            out.write("set ylabel \"Coverage\"\n");
            out.write("set yrange [0:1000]\n");

            // set exon arrows:
            for (Integer border: gene.getLocalExonBoundries()){
                out.write("set arrow from " + border + ", 0 to " +border +",100 lt 5 \n");
            }

            if (gene.getAnnotation().isReverseStrand()){
                out.write("set xrange [] reverse;\n");
            }

            String pngName = dir.getAbsolutePath() + "/" + f.replace(".DoC", ".png");
            //echo "set terminal png;set output 'test.png';plot 'ENST00000084795.DoC' using 4" | gnuplot -persist
            out.write("set output '"+pngName+"';\n");

            // plot command:
            out.write("plot '");
            out.write(dir.getAbsolutePath()+"/"+f);
            out.write("' using 1 notitle with linespoints\n");

            // reset to clear exons for next loop:
            out.write("reset\n");

            // add file to delete list
            tempFiles.add(dir.getAbsolutePath() + "/" + f);
            // break is over max
            count++;
            if (count > MAX_GNUPLOTS) break;
        } // end of genes (and plots)

        out.close();
        gobbler.join();

        for (String file: tempFiles){
            // we can down delete the DoC file for this transcript since we don't need it any more
            File docFile = new File(file);
            docFile.delete();
        }
	}


    /**
     * Reads in the DoC data into an indexed list of coordinates.
     * Then it cycles through the transcripts and look up coverage value for each position in in each transcript.
     * Finally a DoCMetrics obj is created for each transcript and appended to that obj
     *
     * Also, a small DoC file is created for GNUPlot
     *
     * @param docResultFile
     * @param transcripts   after the call to this method, the transcript object contain all the doc results for that transcript
     * @throws IOException
     */
    private void splitDoCResultsByGene(String docResultFile,
            Collection<Transcript> transcripts, String sampleDir, boolean writeSplitsToFile) throws IOException{

        // first load the positions and their values into a hash table
        HashMap<String, String> docMap = indexDoCResults(docResultFile);

        // now we have a map of all the coordinates and their doc result values.
        // now we go through all the transcripts and add the data from the map:
        //ArrayList<TranscriptDoCResult> results = new ArrayList<TranscriptDoCResult>();
        System.out.println("Finding DoC results per transcript");
        Performance perf = new Performance("DoC Results by interval have been mapped back to the transcripts: " + docResultFile);
        int docSplitCount = 0;

        String docMetricsFile = sampleDir+"/transcriptLevelDoCMetrics.txt";
        System.out.println("Writing transcript level DoC metrics to: " + docMetricsFile);
        BufferedWriter out = new BufferedWriter (new FileWriter(docMetricsFile));

        out.write("id\tmeancov\tfiveprime\tthreeprime\t3to5\n");

        int count = 0;
        Performance perfT = new Performance("Processed first 1000 transcripts ",Performance.Resolution.minutes);
        for (Transcript t: transcripts){
            count++;
            if (count % 1000 == 0) {
                System.out.println(perfT);
                perfT = new Performance("Processed " +(count-1000)+" - " + count, Performance.Resolution.minutes);
            }
            TranscriptDoCResult docResult = new TranscriptDoCResult(t.getTranscriptId(),t.getAnnotation().getStrand());
            Collection<int[]> exons = t.exons;
            for (int[] exon: exons) {
                for (int pos = exon[0]; pos <= exon[1]; pos++) {
                    // for each position, we look up the data
                    //construct coordinate from position:

                    String coord = t.getChr() +":" +pos;

                    String dataLine = docMap.get(coord);
                    if (dataLine == null || dataLine.equals("")) {
                        //throw new RuntimeException("No DoC data found for "+t.getTranscriptId()+", " + coord + " in file: " + docResultFile);
                        //System.out.println("Warning, no DoC data found for "+t.getTranscriptId()+", " + coord + " in file: " + docResultFile);
                    } else {
                        docResult.addLine(t.getTranscriptId()+"\t"+t.getAnnotation().getName()+"\t"+dataLine);
                    }
                }
            }
//            System.out.println(perfT);

            t.setDocMetrics(new DoCMetrics(docResult)); // use result obj to create metrics obj

            // this is the last chance to do anything with the actual result values. so we output them for GNUPlot
            if (docSplitCount < MAX_GNUPLOTS){
                docSplitCount++;
                if (writeSplitsToFile) {
                    docResult.writeToFile(sampleDir);
                }
            }
            DoCMetrics mets = t.getDoCMetrics();

            out.write(t.getTranscriptId()); out.write('\t');
            out.write(""+mets.getMeanCoverage()); out.write('\t');
            out.write(""+mets.getFiveEnd400()); out.write('\t');
            out.write(""+mets.getThreeEnd400()); out.write('\t');
            if (mets.getFiveEnd400()>0) {
                float five = mets.getFiveEnd400();
                float three = mets.getThreeEnd400();
                out.write(""+three/five);
            } else {
                out.write("NA");
            }
            out.write('\n');
        }
        out.close();
        System.out.println(perf);
    }


    /**
     * indexes the doc results (coverage) by coordinate
     */
    private HashMap<String, String> indexDoCResults(String docResultFile) throws IOException {
        HashMap<String, String> docMap = new HashMap<String,String>();
        Performance perf = new Performance("Indexing DoC result file: " + docResultFile);
        //System.out.println("(heap is " + Runtime.getRuntime().totalMemory()+")");

        BufferedReader in = new BufferedReader (new FileReader(docResultFile));

        String header = in.readLine();
        String line = in.readLine(); // skipping header
        while (line != null){
            String[] split = line.split("\\t");
            String coord = split[0];
//            String chr = split[0].substring(0, split[0].indexOf(':'));
//            int pos = Integer.parseInt(split[0].substring(split[0].indexOf(':')+1));

            docMap.put(coord,line);
            line = in.readLine();
        }
        in.close();

        System.out.println(perf.toString());
        //System.out.println("(heap is " + Runtime.getRuntime().totalMemory()+")");

        perf = new Performance("Created TranscriptDoCResults for " + docResultFile);
        return docMap;
    }


//    private static Transcript getGeneForPosition(Collection<Transcript> transcripts,String chr, int pos) {
//		Transcript found = null;
//		for (Transcript g: transcripts){
//			if (g.hasPosition(chr,pos)){
//				if (found != null){
//					// found twice
//					System.out.print("Two transcripts have same position");
//					System.out.println("\t"+g.getTranscriptId()+" and " + found.getTranscriptId());
//				}
//				found =  g;
//			}
//		}
//		return found;
//	}

//	/**
//	 * This is the file required by GATK DepthofCoverage which defines the intevals (exons)
//	 * The method just writes one interval per line
//	 *
//	 * @param outfile
//	 * @param intervals
//	 * @throws IOException
//	 */
//	private static void createIntervalList(String outfile,
//			Collection<String> intervals) throws IOException{
//		// output one coordinate per line:
//		BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
//		for (String t: intervals){
//			String coords = t.substring(t.indexOf(' ')+1);
//			StringTokenizer toks =new StringTokenizer(coords,",");
//			while (toks.hasMoreTokens()){
//				out.write(toks.nextToken());
//				out.write('\n');
//			}
//			//System.out.println(coords);
//
//		}
//		out.close();
//	}


	/**
	 * Reads through the expression file to find the MAX-best entries by saving them in a linked list
	 * The file has the form:
     * Gene_Symbol     Ensembl_Transcript_ID   Ensembl_biotype Transcript_Size Total_Transcript_Coverage       Average_Transcript_Coverage     Transcript_Coverage_StdDev
	 *
     * The standard pipline assumes that the first col is id and the 
     * @param infile
     * @param idCol
     *@param expCol @return
     * @throws IOException
     */
    private static TranscriptList readTranscripts(String infile, int cap, int idCol, int expCol) throws IOException{
        BufferedReader in = new BufferedReader(new FileReader (infile));
        LinkedList<Transcript> ll = new LinkedList<Transcript>();

        String line = in.readLine();
        int lineCount=1;
        line = in.readLine(); //skipping header
        while (line!=null){
            String[] split = line.split("\\t");
            double avCov = Double.parseDouble(split[expCol]);

            int i =0;
            if (ll.size()>0){
                float thisAvCov;
                do {
                    thisAvCov =ll.get(i).getExpression();
                    if(thisAvCov < avCov)i++;
                }while (thisAvCov < avCov && i < ll.size());
            }

            ll.add(i,new Transcript(split[idCol],Float.parseFloat(split[expCol])));
            if (ll.size() > cap){
                ll.remove(0);
            }

            line = in.readLine();
            lineCount++;
        }
        in.close();

        // reverse order:
        TranscriptList reverse = new TranscriptList();
        for (int i = ll.size()-1; i >=0; i--){
            reverse.add(ll.get(i));
        }
        return reverse;
    }


    /**
     * all the the lists of transcripts will be found under the path:
     *
     *
     *  OUT_DIR/[low;high;med]expr/sampleId/samplid+med_transcripts.list
     */
    public void createExpressionStratifiedTranscriptLists(int n, String gctFile, float lowerExprCutoff) throws IOException{

        System.out.println("Stratifying Transcripts By Expression");

        GCTFile gct = new GCTFile(gctFile);
        for (RNASeqMetrics.MetricSample b: samples) { // sample id, bam file, notes
            String sampId = b.sampId;

            String[] stratifiedTranscripts = getExprStratifiedTranscripts( n, gct.getTranscripts(sampId), lowerExprCutoff);

            //low expression
            String sampDir = b.getSampleDirectory();
            File dir = new File(b.getExpressionDir(ExpressionLevel.low));
            dir.mkdir();
            String lowExprFile = b.getExpressionList(ExpressionLevel.low);
            IOTools.stringToFile(lowExprFile,stratifiedTranscripts[0]);

            //medium expression
            dir = new File(b.getExpressionDir(ExpressionLevel.medium));
            dir.mkdir();
            String medExprFile = b.getExpressionList(ExpressionLevel.medium);
            IOTools.stringToFile(medExprFile,stratifiedTranscripts[1]);

            //high expression
            dir = new File(b.getExpressionDir(ExpressionLevel.high));
            dir.mkdir();
            String highExprFile =b.getExpressionList(ExpressionLevel.high);
            IOTools.stringToFile(highExprFile,stratifiedTranscripts[2]);

            // //all expression
            // dir = new File(b.getExpressionDir(ExpressionLevel.all));
            // dir.mkdir();
            // String allExprFile =b.getExpressionList(ExpressionLevel.all);
            // IOTools.stringToFile(allExprFile,stratifiedTranscripts[3]);
        }
    }


    /**
     *
     * @return Each String in this array contains a linefeed-separated list of 1000 ids.
     *              [0] = lowly expressed
     *              [1] = expressed around the median
     *              [2] = highly expressed
     *
     * @throws IOException
     */
    private static String[] getExprStratifiedTranscripts(int numberOfTranscripts, ArrayList<GCTFile.GCTEntry> transcripts, float LOWER_EXPR_CUTOFF) throws IOException{
        String[] lists = new String[4];

        ArrayList<GCTFile.GCTEntry> expressed = GCTFile.GCTEntry.sortExpressed(transcripts, LOWER_EXPR_CUTOFF);

        System.out.println("Number of expressed transcripts at this cuttoff: \t"+expressed.size());

        StringBuilder list = new StringBuilder();
        for (int i = 0 ; i < numberOfTranscripts && i < expressed.size(); i ++) {
            list.append(expressed.get(i).getId()).append('\n');
        }
        lists[0] = list.toString();
        list = new StringBuilder();
        int mid = expressed.size()/2;
        int start = Math.max(0, mid-(numberOfTranscripts/2));
        int stop = Math.min(expressed.size(),start+numberOfTranscripts);
        for (int i = start ; i < stop; i ++) {
            list.append(expressed.get(i).getId()).append('\n');
        }
        lists[1] = list.toString();
        list = new StringBuilder();
        for (int i = expressed.size()-1 ; i >= (expressed.size()-numberOfTranscripts) && i >=0; i --) {
            list.append(expressed.get(i).getId()).append('\n');
        }
        lists[2] = list.toString();
        list = new StringBuilder();
        for (int i = 0 ; i < expressed.size(); i ++) {
            list.append(expressed.get(i).getId()).append('\n');
        }

        lists[3] = list.toString();
        return lists;
    }


    public void prepareIntervals(String transcriptGTF, String transcriptTypeField) throws IOException {
        prepareIntervals( ExpressionLevel.low, transcriptGTF, transcriptTypeField);
        prepareIntervals(ExpressionLevel.medium, transcriptGTF,  transcriptTypeField);
        prepareIntervals(ExpressionLevel.high, transcriptGTF,  transcriptTypeField);
        // prepareIntervals(ExpressionLevel.all, transcriptGTF,  transcriptTypeField);  // breaks pre-1.1.9 version
    }

    /**
     * The GATK DepthOfCoverage requiers an interaval list as input. This method creates such interval lists for all
     * samples, but for a given expression level. The GTF is used to determine the intervals. Which transcript to use
     * is determined by the file found in sampleDir/sample_id+exprSuffix todo this way of finding the expression file is nasty
     * can we pack this better into the MetricSample object?
     *
     * @param transcriptGTF
     * @param transcriptTypeField
     * @throws IOException
     */
    public void prepareIntervals(ExpressionLevel level, String transcriptGTF, String transcriptTypeField) throws IOException {

        for (RNASeqMetrics.MetricSample b: samples) { // sample id, bam file, notes

            System.out.println("Expression file for DoC:\t"+ b.getExpressionList(level));
            PerBaseDoC.cleanDir(b.getExpressionDir(level));

            // get MAX best transcripts from expression file
            // e.g. /Volumes/ddeluca/data/DoC/CLL-JP.rna.geneexpression.txt
            System.out.println("Loading transcripts");
            TranscriptList transcripts = Transcript.readTranscriptsFromList(b.getExpressionList(level));

            // use the GTF file to load the intervals, as well as some annotation information
            Performance perf = new Performance("Interval Loading");
            System.out.println("Preparing intervals for " + transcripts.size() + " transcripts");

            TranscriptList fullTrans = TranscriptList.loadGTF(transcriptGTF,transcriptTypeField);

            transcripts.loadIntervalsAndAnnotation(fullTrans);
            System.out.println(perf);

            // create interval list for GATK
            System.out.println("Creating interval list");

            Transcript.createIntervalFile(transcripts,b.getIntervalFile(level));
        }
    }


    public void runDoC(String refGenomeFile) throws Exception {
        runDoC(refGenomeFile,null,null,false); // no downsampling
    }


    public void runDoC(String refGenomeFile, String dSamplingTarget, String additionalFlags, boolean runOnAll) throws Exception {
        runDoC(refGenomeFile, ExpressionLevel.low, dSamplingTarget, additionalFlags);
        runDoC(refGenomeFile,ExpressionLevel.medium, dSamplingTarget, additionalFlags);
        runDoC(refGenomeFile,ExpressionLevel.high, dSamplingTarget, additionalFlags);
        // if (runOnAll) {  // breaks pre-1.1.9 version
        //     runDoC(refGenomeFile,ExpressionLevel.all, dSamplingTarget, additionalFlags);
        // }
    }


    public void runDoC(String refGenomeFile,ExpressionLevel level, String dSamplingTarget, String additionalFlags_) throws Exception {

        for (RNASeqMetrics.MetricSample b: samples){ // sample id, bam file, notes

            String intervalFile = b.getIntervalFile(level);
            // run DoC
            String docResults = b.getDoCResultsFile(level);
            String finalFlags = additionalFlags_;
            //-rf Downsampling -numReads 2000 -targetReads 1000
            String dSampleOption = null;
            if (dSamplingTarget !=null){
                int totReads = b.getTotalReads();
                System.out.println("Coverage Metrics are being downsampled from " +totReads +" + for sample: " + b.sampId);
                dSampleOption = "-rf Downsampling -numReads " + totReads + " -targetReads " + dSamplingTarget;
                if (finalFlags == null) {
                    finalFlags= dSampleOption;
                } else {
                    finalFlags += " " + dSampleOption;
                }
            }

            if (!noGatk) {
                GATKTools.runDoC(refGenomeFile, b.getBamFileOrList(), docResults, intervalFile, false, finalFlags);
            } else {
                System.out.println("GATK DoC suppressed");
            }
        }
    }


    /**
     *
     * @param detailedHTML
     * @param gcFile         can be null, in which case it is ignored
     * @param endLength
     * @throws IOException
     */
    public void createReports(String transcriptGTF, String transcriptTypeField,boolean detailedHTML,String gcFile, String endLength) throws IOException {
        createReports(transcriptGTF, transcriptTypeField, ExpressionLevel.low, detailedHTML, gcFile, endLength);
        createReports(transcriptGTF, transcriptTypeField, ExpressionLevel.medium, detailedHTML, gcFile, endLength);
        createReports(transcriptGTF, transcriptTypeField, ExpressionLevel.high, detailedHTML, gcFile, endLength);
        // createReports(transcriptGTF, transcriptTypeField, ExpressionLevel.all, detailedHTML, gcFile, endLength);
    }

    private void createReports(String transcriptGTF, String transcriptTypeField,ExpressionLevel level,  boolean detailedHTML, String gcFile, String endLength) throws IOException {

        interpretAndOutputDoCResults(transcriptGTF, transcriptTypeField, level, detailedHTML, gcFile, endLength);

        // creates a single file per expression level with the data used in the coverage plots
        agglomerateCoveragePlotData(level);
        agglomerateCoveragePlotNormalizedData(level);

        // creates a single file per expression level with the data used in the gap length histograms
        agglomerateGapLengthHisograms(level);

        // creates an HTML page for each sample and coverage stratification
        createDoCPages(level, detailedHTML);
    }

    /**
     * this is where the real analysis of the DoC output occus 
     *
     * @param transcriptGTF
     * @param transcriptTypeField
     * @param level
     * @param detailedHTML
     * @param gcFile
     * @param endLength
     * @throws IOException
     */
    private void interpretAndOutputDoCResults(String transcriptGTF, String transcriptTypeField, ExpressionLevel level, boolean detailedHTML, String gcFile, String endLength) throws IOException {
        TranscriptList fullTrans = TranscriptList.loadGTF(transcriptGTF,transcriptTypeField);
        for (RNASeqMetrics.MetricSample b: samples){ // sample id, bam file, notes

            System.out.println("Loading transcripts");
            TranscriptList transcripts = Transcript.readTranscriptsFromList(b.getExpressionList(level));

            transcripts.loadIntervalsAndAnnotation(fullTrans);
            if (gcFile !=null) {
                //populate transcripts with GC content information
                transcripts.populateGC(gcFile);
            }

            //read DoC results and split them into files by gene (for gnuplot)
            System.out.println("Splitting intervals into transcript-oriented DoC files");


            Performance perf = new Performance("Mapped intervals back to transcripts",Performance.Resolution.seconds);
            splitDoCResultsByGene(b.getDoCResultsFile(level),transcripts, b.getExpressionDir(level), detailedHTML); // outDir is just for making GNUPlot input files
            System.out.println(perf);

            //create data on the averaage coverage at each position relative to 3' end
            calculateAverageCoverageByPosition(transcripts, b.getMeanCoverageByPosFile(level));

            calculateAverageCoverageByNormalizedPosition(transcripts, b.getMeanCoverageByPosNormalizedFile(level));
            calculateGapLengthHistogram(transcripts,b.getGapLengthHistogramFile(level));

            // at this point the transcripts contain the doc Result info
            transcripts.toFile(b.getDoCTranscriptsResultsFile(level), endLength, gcFile != null, detailedHTML);

            if (detailedHTML) {
                //create GNU Plots:
                System.out.println("Creating Plots");
                try {
                    createGNUPlots(transcripts, b.getExpressionDir(level));
                } catch(IOException e){
                    System.out.println("\nCould not run GNUPlot. Coverage plots for individual transcripts will not be created.\n");
                } catch (RuntimeException e){
                    System.out.println("\nCould not run GNUPlot. Coverage plots for individual transcripts will not be created.\n");
                } catch (InterruptedException e) {
                    System.out.println("\nGNUPlot Interrupted. Coverage plots for individual transcripts will not be created.\n");
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
                // create html report:
                System.out.println("Creating HTML Report");
                createDetailedHTML(b.getExpressionDir(level), transcripts, gcFile!=null);
            }

            transcripts.toSummaryFile(b.getDoCTranscriptsSummaryResultsFile(level));
        }
    }


    /**
     * the output file will look like OUT_DIR+"/lowexpr/meanCoverage.txt"
     *
     * @throws IOException
     */
    private void agglomerateCoveragePlotData(ExpressionLevel level) throws IOException {

        //set up data for plot
        XYSeriesCollection xySeriesCollection = new XYSeriesCollection();

        ArrayList<String>[] table= new ArrayList[samples.size()];

        BufferedWriter out = new BufferedWriter(new FileWriter(this.outDir + "/meanCoverage_"+level+".txt" ));
        int length = 0;
        for (int i = 0 ; i < table.length; i++){
            String sampId = samples.get(i).sampId;
            table[i] = IOTools.fileToList(samples.get(i).getMeanCoverageByPosFile(level), null);
            length = table[i].size();
            out.write(sampId);
            out.write('\t');

            //add a new series to the plot
            XYSeries xySeries = new XYSeries(sampId);
            for(int r =0;r<table[i].size();r++) {
                double value = Double.parseDouble(table[i].get(r));
                xySeries.add(r+1, value);
            }

            xySeriesCollection.addSeries(xySeries);
        }
        out.write('\n');

        for (int row = 0 ; row < length; row++) {
            for (int sampDex = 0 ; sampDex < table.length; sampDex++){
                out.write(table[sampDex].get(row)); out.write('\t');
            }
            out.write('\n');
        }
        out.close();

        //save mean coverage plot
        System.setProperty("java.awt.headless", "true");
        JFreeChart meanCoveragePlot = ChartFactory.createXYLineChart
                              (null, "Distance from 3'", "Mean Coverage", xySeriesCollection, PlotOrientation.VERTICAL, true, true, false);
        saveImage(meanCoveragePlot, new File(this.outDir + "/meanCoverage_"+level+".png"));
    }


    /**
     * the output file will look like OUT_DIR+"/lowexpr/meanCoverage.txt"
     *
     * @throws IOException
     */
    private void agglomerateCoveragePlotNormalizedData(ExpressionLevel level) throws IOException {

        //set up data for plot
        XYSeriesCollection xySeriesCollection = new XYSeriesCollection();

        ArrayList<String>[] table= new ArrayList[samples.size()];

        BufferedWriter out = new BufferedWriter(new FileWriter(this.outDir + "/meanCoverageNorm_"+level+".txt" ));
        int length = 0;
        for (int i = 0 ; i < table.length; i++){
            String sampId = samples.get(i).sampId;
            table[i] = IOTools.fileToList(samples.get(i).getMeanCoverageByPosNormalizedFile(level), null);
            length = table[i].size();
            out.write(sampId);
            out.write('\t');

            //add a new series to the plot
            XYSeries xySeries = new XYSeries(sampId);
            for(int r =0;r<table[i].size();r++) {
                double value = Double.parseDouble(table[i].get(r));
                xySeries.add(r+1, value);
            }

            xySeriesCollection.addSeries(xySeries);
        }
        out.write('\n');

        for (int row = 0 ; row < length; row++) {
            for (int sampDex = 0 ; sampDex < table.length; sampDex++){
                out.write(table[sampDex].get(row)); out.write('\t');
            }
            out.write('\n');
        }
        out.close();

        //save mean coverage plot
        System.setProperty("java.awt.headless", "true");
        JFreeChart meanCoveragePlot = ChartFactory.createXYLineChart
                (null, "Percentage of Transcript Length (5' to 3')", "Mean Coverage", xySeriesCollection, PlotOrientation.VERTICAL, true, true, false);
        saveImage(meanCoveragePlot, new File(this.outDir + "/meanCoverageNorm_"+level+".png"));
    }


    /**
     * the output file will look like OUT_DIR+"/lowexpr/gapLengthHist.txt"
     *
     * @throws IOException
     */
    private void agglomerateGapLengthHisograms(ExpressionLevel level) throws IOException {

        //set up data for plot
        XYSeriesCollection xySeriesCollection = new XYSeriesCollection();
        ArrayList<String>[] table= new ArrayList[samples.size()];

        BufferedWriter out = new BufferedWriter(new FileWriter(this.outDir + "/gapLengthHist_"+level+".txt" ));
        int length = 0;
        for (int i = 0 ; i < table.length; i++) {
            String sampId = samples.get(i).sampId;
            table[i] = IOTools.fileToList(samples.get(i).getGapLengthHistogramFile(level), null);
            length = table[i].size();
            out.write(sampId);
            out.write('\t');

            //add a new series to the plot
            XYSeries xySeries = new XYSeries(sampId);
            for(int r =0;r<table[i].size();r++) {
                double value = Double.parseDouble(table[i].get(r));
                xySeries.add(r+1, value);
            }

            xySeriesCollection.addSeries(xySeries);
        }
        out.write('\n');

        for (int row = 0 ; row < length; row++) {
            for (int sampDex = 0 ; sampDex < table.length; sampDex++) {
                out.write(table[sampDex].get(row)); out.write('\t');
            }
            out.write('\n');
        }
        out.close();

        //save mean coverage plot
        System.setProperty("java.awt.headless", "true");
        JFreeChart plot = ChartFactory.createXYLineChart
                (null, "Gap Length'", "Number of Gaps", xySeriesCollection, PlotOrientation.VERTICAL, true, true, false);
        saveImage(plot, new File(this.outDir + "/gapLengthHist_"+level+".png"));
    }


    public static void saveImage(JFreeChart chart, File file) {

        ChartPanel panel = new ChartPanel(chart);
        Dimension size = panel.getPreferredSize();

        try {
            ChartUtilities.saveChartAsPNG(file, chart, size.width, size.height);
        }
        catch(IOException io) {
            io.printStackTrace();
            System.err.println("Unable to save the image file " + file);
        }
    }

    /**
     * This creates an index.html for each sample (lane)
     *
     * The table is a list of the transcripts for this lane and the transcript-specific stats
     *
     * This file will be located in each of the highexpr, lowexpr and medexpr folders
     * @throws IOException
     */
    private void createDoCPages(ExpressionLevel level, boolean link) throws IOException {

        for (RNASeqMetrics.MetricSample samp: this.samples) {
            String sampleId = samp.sampId;

            // output all the transcripts in a big table
            BufferedWriter out = new BufferedWriter(new FileWriter(samp.getExpressionDir(level)+"/index.html"));
            out.write("<html>\n<head><title>"+sampleId+"</title>\n\t"+RNASeqMetrics.getStyle()+"\n\t</head>\n<body>");
            out.write("<h1>"+sampleId+"</h1\n");
            out.write("<p>\n<h2>Data Files</h2>\n");
            out.write("<dl>\n");
            if((new File(samp.getExpressionDir(level)+"../" + sampleId+".metrics.txt")).exists()) {
                out.write("     <dt><a href='../../"+sampleId+"/"+sampleId+".metrics.txt'>Read Count Metrics</a></dt>");
                out.write("     <dd>Tab delimited text file of the read count-based metrics for this sample.</dd>\n");
            }

            System.out.println("Library size link path" + new File(samp.getExpressionDir(level) + "../"));
            if((new File(samp.getExpressionDir(level)+"../" + sampleId+".libraryComplexity.txt")).exists()) {
                out.write("     <dt><a href='../../"+sampleId+"/"+sampleId+".libraryComplexity.txt'>Library Complexity Metrics</a></dt>");
                out.write("     <dd>Contains estimated library size, number of read pairs, unpaired reads, etc.</dd>\n");
            }

            System.out.println("Absolute path to metrics.tmp.txt: " + (new File("../" + samp.getExpressionDir(level)+"/" + sampleId+".metrics.tmp.txt.rpkm.gct").getAbsolutePath()));
            if ((new File(samp.getExpressionDir(level)+"../" + sampleId+".metrics.tmp.txt.rpkm.gct")).exists()) {
                out.write("     <dt><a href='../../"+sampleId+"/"+sampleId+".metrics.tmp.txt.rpkm.gct'>Expression File</a></dt>");
                out.write("     <dd>GCT file containing the expression levels of this sample in RPKM</dd>\n");
            }
            out.write("     <dt><a href='perBaseDoC.out'>Coverage per Base</a></dt>");
            out.write("     <dd>GATK Depth Of Coverage output for each genomic position</dd>\n");
            out.write("     <dt><a href='perBaseDoC.out.sample_summary'>GATK Depth of Coverage Summary</a></dt>");
            out.write("     <dd>Top level metrics file created by GATK after running Depth of Coverage</dd>\n");

            out.write("</dl>\n");
            out.write("</p>\n");

            if (link) {
                out.write("<p>");
                out.write("<a href='./plots.html'>View All Plots</a>");
                out.write("</p>");
            }

            out.write("<h2>Transcript Details</h2>\n");

            out.write(IOTools.tableFileToHTML(samp.getDoCTranscriptsResultsFile(level)));
            out.write("<p class='tLegend'>\n");
            out.write("<b>CV</b> is coefficient of variation (standard deviation of coverage divided by mean coverage). ");
            out.write("<b>5'</b> and <b>3'</b> values are per-base coverage. <b>Norm</b> denotes that the end coverage is divided by the mean coverage for that transcript.");
            out.write(" <b>Gaps</b> are defined as a region of 5 or more bases having zero coverage. ");
            out.write(" <b>Gap %</b> is the total cumulative gap length divided by the total cumulative transcript lengths.");
            out.write("\n</p>\n");

            out.write("<h5>"+new Date()+"</h5>\n");
            out.write("</body>\n</html>");
            out.close();
        }
    }

}
