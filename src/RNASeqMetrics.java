package org.broadinstitute.cga.rnaseq;

import org.apache.commons.cli.*;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math.stat.correlation.SpearmansCorrelation;
import org.broadinstitute.cga.tools.*;

import java.io.*;
import java.util.*;


public class RNASeqMetrics {
    public static String VERSION = "v1.1.9 06/26/16";

    static String USAGE = "java -jar RNA-SeQC.jar <args>";
    private static  float LOWER_GC_CUTOFF = 0.375f;
    private static  float UPPER_GC_CUTOFF = 1f-LOWER_GC_CUTOFF;
    public static final float LOWER_EXPR_CUTOFF =  1f;

    private boolean noMetrics;
    int MAX = 1000;
    String SAMPLE_FILE;
    String TRANSCRIPT_MODEL;
    String REF_GENOME;
    String OUT_DIR;
    boolean details;
    String gatkFlags = null;
    String transcriptTypeField ;
    String endLength;
    String refGeneFile;
    String rRNAIntervals;
    String rRNAIntervalsCreationLocation;
    String rRNAFile;
    private String theseSamplesGCT =null;
    private String dSampleTarget;
    private String gcFile;
    private boolean doStrat;
    private String referenceGCTForCorr;
    ArrayList<MetricSample> bams;
    private boolean createRefGeneFromGTF = false;
    private boolean filterGCT;
    private String bwa;
    private boolean singleEnd;
    private boolean isStrat;
    private boolean docAll = false;
    private int rRNAdSampleTarget = 1000000;
    private boolean gapLengthDistribution = false;

    ArrayList<HashMap<String,String>> metricTracker = new ArrayList<HashMap<String,String>>();
    private boolean skiprRNA = false;
    private boolean runDoC = true;
    private boolean strictMode = false;


    private static Options setupCliOptions() {
        Options opts = new Options();
        // first the non-argument, optional parameters:
        opts.addOption("noDoC", false, "Suppresses GATK Depth of Coverage calculations.");
        opts.addOption("noReadCounting", false,  "Suppresses read count-based metrics.");
        opts.addOption("transcriptDetails", false,  "Provide an HTML report for each transcript.");
        opts.addOption("singleEnd", false,  "This BAM contains single end reads.");
        opts.addOption("gld", false,  "Gap Length Distribution: plots the distribution of gap length.");
        opts.addOption("strictMode",false, "Applies strict filtering before quantifying transcripts");
        opts.addOption("fullDoC",false, "Creates depth-based metrics on every expressed gene");

        // required parameters:
        Option [] required = { 
                new Option("s", "Sample File (text-delimited description of samples and their bams)."),
                new Option("t", "GTF File defining transcripts (must end in '.gtf')."),
                new Option("r", "Reference Genome in fasta format."),
                new Option("o", "Output directory (will be created if doesn't exist).")
        };
        // set all other options are required, and having one argument and add to Options object
        for (Option o:required) {
            o.setRequired(true);
            o.setArgs(1);
            opts.addOption(o);
        }
        // optional parameters with one agument:
        Option [] optional = {  new Option("n", "Number of top transcripts to use."),
                new Option("d", "Perform downsampling to the given number of reads."),
                new Option("e", "Change the definition of a transcripts end (5' or 3') to the given length. " +
                        "           (10, 50, 100 are acceptable values)"),
                new Option("rRNA", "intervalFIle for rRNA loci (must end in .list)"),
                new Option("corr", "GCT file for expression correlation comparison"),
                new Option("gc", "File of transcript id <tab> gc content. Used for stratification."),
                new Option("strat", "Stratification options: current supported option is 'gc'"),
                new Option("expr", "Uses provided GCT file for expression values instead of on-the-fly RPKM calculation"),
                new Option("BWArRNA", false,  "Use an on the fly BWA alignment for estimating rRNA content. The value should be the rRNA reference fasta."),
                new Option("ttype", "The column in gtf to use to look for rRNA transcript type"),
                new Option("bwa", "Path to BWA, which should be set if it's not in your path and BWArRNA is used."),
                new Option("rRNAdSampleTarget", "Downsamples to calculate rRNA rate more efficiently. Default is 1 million. Set to 0 to disable."),
                new Option("gcMargin", "Used in conjunction with '-strat gc' to specify the percent gc content to use as boundaries. E.g. .25 would set a lower cutoff of 25% and an upper cutoff of 75%."),
                new Option("gatkFlags", "A string of flags that will be passed on to the GATK"),
        };

        for (Option o:optional) {
            o.setRequired(false);
            o.setArgs(1);
            opts.addOption(o);
        }
        return opts;
    }


    private static void printHelp(Options opt) {
        HelpFormatter f = new HelpFormatter();
        f.printHelp(USAGE, opt);
    }


    private static void checkArgs(CommandLine cl) throws Exception{
        if (cl.hasOption("e")){
            String l = cl.getOptionValue("e");
            // end length
            if (!l.equals("200") && !l.equals("50") && !l.equals("100")){
                throw new Exception("Acceptable values for e are 50, 100 or 200");
            }
        }
        if(!cl.getOptionValue("t").toLowerCase().endsWith(".gtf")){
            throw new Exception("Transcript file must be in GTF format (and end with the '.gtf' extension); file provided: " + cl.getOptionValue("t"));
        }
    }

    /**
     *
     * @param argz
     */
    public static void main(String[] argz) {
        Performance perf = new Performance("RNA-SeQC Total Runtime", Performance.Resolution.minutes);
        Options opts = setupCliOptions();
        try {
            execute(argz);
            System.out.println("Finished Successfully.");
            System.out.println(perf);
        } catch (MissingOptionException moe) {
            System.err.println(moe.getMessage());
            printHelp(opts);
            System.exit(1);
        } catch (ParseException moe) {
            moe.printStackTrace();
            printHelp(opts);
            System.exit(2);
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println(perf);
            System.exit(3);
        }
    }


    private static void execute(String[] argz) throws Exception{
        System.out.println("RNA-SeQC " + VERSION);
        Options opts = setupCliOptions();
        CommandLineParser parser = new PosixParser();
        CommandLine cl = parser.parse(opts, argz);
        checkArgs(cl);

        RNASeqMetrics metrics = new RNASeqMetrics();
        metrics.prepareOptionsFromCommandLine(cl);
        metrics.prepareFiles(); // makes the rout directory, does the GTF to refgene conversion, creates rRNA intervals, filters user defined GCT to match GTF
        metrics.runMetrics();

        if (metrics.runDoC && cl.hasOption("strat") && !cl.getOptionValue("strat").equals("none")) {
            stratification(cl, argz);
        }
    }


    private void prepareFiles() throws IOException {

        // create directory
        File dir = new File(OUT_DIR);
        dir.mkdir();

        // load samples (sample id, expression file, bam file, notes)
        bams = MetricSample.readInSamples(SAMPLE_FILE, OUT_DIR);
        if (bams.size() == 0){
            throw new RuntimeException("The list of samples was empty. There was something wrong with your -s tag.");
        }

        // create refGene file
        try {
            createRefGeneAndRRNAFiles(TRANSCRIPT_MODEL, REF_GENOME, refGeneFile, rRNAIntervalsCreationLocation, transcriptTypeField);
        } catch (RuntimeException e) {
            e.printStackTrace();
            System.out.println("No information for rRNA available. Continuing without rRNA calculations. (Using the -BWArRNA flag for best results)");
            skiprRNA = true;
        }
        if (filterGCT) {
            // this expressino table might have values rows that do not correspond to the transcript def.
            // they must be filtered out
            System.out.println("Filtering GTF file to correspond to GCT file.");
            String filteredGCT = OUT_DIR+"/filtered_expression.gtf";
            createFilteredGCT(TRANSCRIPT_MODEL,theseSamplesGCT, filteredGCT);
            theseSamplesGCT = filteredGCT;
        }
        
        // initialize metric tracker for this number of samples:
        for (int i = 0 ; i < bams.size(); i++) {
            metricTracker.add(new HashMap<String, String>());
        }
    }


    private void runMetrics() throws Exception {

        String metricTable = null;
        if (!noMetrics) {
            ReadCountMetrics rcMetrics = new ReadCountMetrics(bams,OUT_DIR, LOWER_EXPR_CUTOFF);
            rcMetrics.runReadCountMetrics(REF_GENOME,refGeneFile, rRNAFile, this.singleEnd, this.bwa, this.rRNAdSampleTarget, this.skiprRNA, this.gatkFlags, this.strictMode);

            //merge GCT files into one
            String gtc = rcMetrics.getCombinedGCTFileName();
            if (theseSamplesGCT == null) theseSamplesGCT = gtc;
            // load and ouput metrics:
            metricTable = rcMetrics.getMetricsHTML(metricTracker);
        } else {
            System.out.println("Metrics suppressed");
        }

        String corrTable = null;
        if (referenceGCTForCorr !=null || bams.size() > 1) {
            // run R to find the correlation between these rpkms and a given (e.g. affy) microarray
            String[] corFile = correlationComparison(referenceGCTForCorr, theseSamplesGCT, OUT_DIR);

            corrTable = "\n<HR WIDTH=\"100%\" COLOR=\"#9090900\"  SIZE=\"3\">\n";
            corrTable+= ("\n<h2>Correlation Analysis</h2>\n");

            if (referenceGCTForCorr !=null){
                corrTable+= createCorrHTMLTable(bams,corFile[0]); // turn corrTable into html
            }
            if (bams.size() > 1){
                corrTable += "\n<h3>Spearman Correlation Matrix</h3>\n";
                corrTable +=createCorrHTMLTableMatrix(bams,corFile[1]);
                corrTable += "\n<h3>Pearson Correlation Matrix</h3>\n";
                corrTable +=createCorrHTMLTableMatrix(bams,corFile[2]);
            }
        }

        // create individual depth of coverage reports per sample
        if (runDoC) {
            PerBaseDoC doc = new PerBaseDoC(bams, OUT_DIR);

            // this step slices up the transcripts by expression level
            doc.createExpressionStratifiedTranscriptLists(MAX, theseSamplesGCT,LOWER_EXPR_CUTOFF);
            doc.prepareIntervals(TRANSCRIPT_MODEL,transcriptTypeField);
            doc.runDoC(REF_GENOME, dSampleTarget, gatkFlags, docAll);

            // Create Index.html
            doc.createReports(TRANSCRIPT_MODEL,transcriptTypeField, details, this.gcFile, endLength);
            createTopLevelIndex(bams,MAX, metricTable, corrTable,endLength, OUT_DIR+"/index.html");
            if (!isStrat) createTopLevelIndex(bams,MAX, metricTable, corrTable,endLength, OUT_DIR+"/report.html");
        }
        if (isStrat) {
            createTSV(bams, OUT_DIR+"/metrics_strat.tsv");
        } else {
            createTSV(bams, OUT_DIR+"/metrics.tsv");
        }
    }


    private void createTSV(ArrayList<MetricSample> bams, String outfile) throws IOException{
        BufferedWriter out = new BufferedWriter(new FileWriter(outfile));

        // headers
        Set<String> fields = metricTracker.get(0).keySet();
        out.write("Sample\tNote");
        for (String f: fields){
            out.write('\t');
            out.write(f);
        }
        out.write("\n");

        //content
        for (int i = 0 ; i < bams.size(); i++) {
            MetricSample samp = bams.get(i);
            HashMap<String,String> metrics = metricTracker.get(i);
            out.write(samp.sampId);
            out.write("\t");
            out.write(samp.notes);
            for (String f: fields) {
                out.write("\t");
                out.write(metrics.get(f));
            }
            out.write("\n");
        }
        out.close();
    }


    private void prepareOptionsFromCommandLine(CommandLine cl) throws Exception {
        // GET ARGUMENTS****************************
        runDoC = !cl.hasOption("noDoC");
        noMetrics = cl.hasOption("noReadCounting");
        if (cl.hasOption("n")) {
            MAX = Integer.parseInt(cl.getOptionValue("n"));
        }
        SAMPLE_FILE = cl.getOptionValue("s");
        TRANSCRIPT_MODEL = cl.getOptionValue("t");
        REF_GENOME = cl.getOptionValue("r");
        OUT_DIR = cl.getOptionValue("o");
        details = cl.hasOption("transcriptDetails");
        transcriptTypeField = cl.getOptionValue("ttype");
        endLength = "200";
        if (cl.hasOption("e")) endLength = cl.getOptionValue("e");
        if (cl.hasOption("rRNAdSampleTarget")) rRNAdSampleTarget = Integer.valueOf(cl.getOptionValue("rRNAdSampleTarget"));
        if (cl.hasOption("gcMargin")) {
            LOWER_GC_CUTOFF = Float.valueOf(cl.getOptionValue("gcMargin"));
            UPPER_GC_CUTOFF = 1f-LOWER_GC_CUTOFF;
        } 
        if (cl.hasOption("gatkFlags")) {
            this.gatkFlags = cl.getOptionValue("gatkFlags");
            System.out.println("Additional GATK flags provided: " + this.gatkFlags);
        }
        this.strictMode = cl.hasOption("strictMode");
        this.docAll = cl.hasOption("fullDoC");
        // ******************************************

        // create directory
        File dir = new File(OUT_DIR);
        dir.mkdir();

        rRNAIntervals= cl.getOptionValue("rRNA");
        rRNAIntervalsCreationLocation = null;
        if (TRANSCRIPT_MODEL.toLowerCase().endsWith(".gtf")) {
            this.createRefGeneFromGTF = true;
            refGeneFile = OUT_DIR +"/refGene.txt";

            if (rRNAIntervals == null && !cl.hasOption("BWArRNA")) {
                System.out.println("Creating rRNA Interval List based on given GTF annotations");
                rRNAIntervalsCreationLocation = OUT_DIR +"/rRNA_intervals.list";
                rRNAIntervals = rRNAIntervalsCreationLocation;
            }
        } else {
            refGeneFile = TRANSCRIPT_MODEL;
            //todo if refGene is provided directly this implementation requires that the rRNAIntervals also
            // be provided, because we are not extracting them from the provided refgene file.
        }

        if (cl.hasOption("strat") && cl.getOptionValue("strat").equals("none") && !cl.hasOption("corr")){
            System.out.println("Suppressing Read Count Metrics within Recursive Call.");
            noMetrics = true; // if there's no correlation, and we're in a recursive call, then we don't need to calculate these metrics.
        }

        if( cl.hasOption("BWArRNA")) {
            rRNAFile = cl.getOptionValue("BWArRNA");
        } else {
            rRNAFile = rRNAIntervals;
        }

        filterGCT = false;
        if (cl.hasOption("expr")) {
            filterGCT = true;
            // this expressino table might have values rows that do not correspond to the transcript def.
            // they must be filtered out
            theseSamplesGCT = cl.getOptionValue("expr");
        }

        referenceGCTForCorr =  cl.getOptionValue("corr");

        dSampleTarget = null;
        if (cl.hasOption("d") && !cl.hasOption("noReadCounting")) {
            dSampleTarget = cl.getOptionValue("d");
        }

        this.gcFile = cl.getOptionValue("gc");
        this.doStrat = cl.hasOption("strat") && !cl.getOptionValue("strat").equals("none");
        this.isStrat = cl.hasOption("strat") && cl.getOptionValue("strat").equals("none");
        this.singleEnd = cl.hasOption("singleEnd");
        this.bwa  = cl.getOptionValue("bwa");
        this.gapLengthDistribution = cl.hasOption("gld");
    }


    private static void createFilteredGCT(String transcriptModel, String expr, String filteredFileName) throws IOException {
        if (!transcriptModel.endsWith("gtf")) throw new RuntimeException("transcript model must be gtf");

        Set<String> modelIds = Transcript.loadGTFIds(transcriptModel)  ;

        GCTFile gct =new GCTFile(expr);
        GCTFile filtered = gct.filterById(modelIds);
        filtered.toFile(filteredFileName);
    }


    private static void stratification(CommandLine cl, String[] argz) throws Exception{
        String TRANSCRIPT_MODEL = cl.getOptionValue("t");
        String OUT_DIR = cl.getOptionValue("o");
        // by GC content
        if (cl.getOptionValue("strat").contains("gc")){  // stratification only supported when GTF file (not refgene) is provided

            HashSet<String>[] stratifiedTranscripts = getGCStratifiedTranscripts(cl.getOptionValue("gc"),2f);

            System.out.println("Copying transcript model in to GC stratifications");
            File dir = new File(OUT_DIR+"/gc");
            dir.mkdir();

            //low GC
            String lowGCFile = OUT_DIR+"/gc/lowgc.gtf";
            Transcript.copyFiltered(TRANSCRIPT_MODEL,lowGCFile, stratifiedTranscripts[0]);
            //low GC
            String midGCFile = OUT_DIR+"/gc/medgc.gtf";
            Transcript.copyFiltered(TRANSCRIPT_MODEL,midGCFile, stratifiedTranscripts[1]);

            //low GC
            String highGCFile = OUT_DIR+"/gc/highgc.gtf";
            Transcript.copyFiltered(TRANSCRIPT_MODEL,highGCFile, stratifiedTranscripts[2]);

            // call this main method for each one:
            mainRecursive(lowGCFile,OUT_DIR+"/gc/low",argz);
            mainRecursive(midGCFile,OUT_DIR+"/gc/mid",argz);
            mainRecursive(highGCFile,OUT_DIR+"/gc/high",argz);
        }
    }


    private static void mainRecursive(String gtfFile, String outdir, String[] argz) throws Exception{
        String[] newArgs = new String[argz.length];
        for (int i = 0 ; i <newArgs.length; i++){
            newArgs[i] = argz[i];
        }

        String oldOutDir = null;
        boolean userRRNA = false;
        boolean userExpr = false;
        for (int i= 0 ; i < newArgs.length; i++){
            if (newArgs[i].equals("-t")){
                newArgs[i+1] = gtfFile;

            }
            if(newArgs[i].equals("-o")){
                oldOutDir = newArgs[i+1];
                newArgs[i+1] = outdir;
            }
            if(newArgs[i].equals("-strat")){
                newArgs[i+1] = "none"; // prevents endless recursion
            }
            if(newArgs[i].equals("-rRNA")){
                userRRNA = true;
            }
            if(newArgs[i].equals("-expr")){
                userExpr = true;
            }
        }

        ArrayList<String> additionalArgs = new ArrayList<String>();
        additionalArgs.add("-noReadCounting");
        if (!userRRNA){
            additionalArgs.add( "-rRNA");
            additionalArgs.add( oldOutDir+"/rRNA_intervals.list");
        }
        if (!userExpr){
            additionalArgs.add("-expr");
            additionalArgs.add(oldOutDir+"/genes.rpkm.gct");
        }
        String[] newArgs2 = new String[newArgs.length+additionalArgs.size()];
        int i = 0;
        for ( ; i < newArgs.length; i++){
            newArgs2[i] = newArgs[i];
        }
        for (String a: additionalArgs){
            newArgs2[i] = a;
            i++;
        }

        System.out.println("Stratifying transcripts with file: \t"+gtfFile);
        System.out.println("\tArguments in this stratification \t" + Arrays.toString(newArgs2));
        RNASeqMetrics.execute(newArgs2);
    }


    private static HashSet<String>[] getGCStratifiedTranscripts(String gcFileName, Float sigmas) throws IOException {
        // read in the ids and GC contents
        BufferedReader in = new BufferedReader (new FileReader(gcFileName));

        String line = in.readLine(); // no header
        ArrayList<String> ids = new ArrayList<String>();
        MathList gcVals = new MathList();
        while (line != null) {
            String [] split = line.split("\\t");
            ids.add(split[0]);
            gcVals.add(Float.valueOf(split[1]));
            line = in.readLine();
        }

        // find low GC transcripts
        HashSet<String> lowGCIds = new HashSet<String>();

        float lowerCutoff = LOWER_GC_CUTOFF; //mean - (sigma * sigmas);
        System.out.println("Lower Bound Z score cutoff: " + lowerCutoff);
        for (int i = 0 ; i < gcVals.size(); i++){
            if (gcVals.get(i) < lowerCutoff){
                lowGCIds.add(ids.get(i));
            }
        }

        System.out.println("\tPercentile:\t" + (float)lowGCIds.size()/(float)gcVals.size());
        System.out.println("\tTotal Transcripts in this percentile: " + lowGCIds.size());

        if (lowGCIds.size() == 0){
            throw new RuntimeException("No Transcripts for lower GC cuttoff: try increasing the value of -gcMargin.");
        }

        // find high GC transcripts
        HashSet<String> highGCIds = new HashSet<String>();

        float upperCutoff = UPPER_GC_CUTOFF ; // mean + (sigma * sigmas);
        System.out.println("Upper Bound Z score cutoff: " + upperCutoff);
        for (int i = 0 ; i < gcVals.size(); i++){
            if (gcVals.get(i) > upperCutoff){
                highGCIds.add(ids.get(i));
            }
        }
        System.out.println("\tPercentile:\t" + (float)highGCIds.size()/(float)gcVals.size());
        System.out.println("\tTotal Transcripts in this percentile: " + highGCIds.size());
        if (highGCIds.size() == 0){
            throw new RuntimeException("No Transcripts for upper GC cuttoff: try increasing the value of -gcMargin.");
        }

        // find high GC transcripts
        HashSet<String> midGCIds = new HashSet<String>();
        for (int i = 0 ; i < gcVals.size(); i++){
            float gc = gcVals.get(i);
            if (gc >= lowerCutoff && gc <= upperCutoff){
                midGCIds.add(ids.get(i));
            }
        }
        System.out.println("\tMiddle Percentile:\t" + (float)midGCIds.size()/(float)gcVals.size());
        System.out.println("\tTotal Transcripts in this percentile: " + midGCIds.size());

        HashSet<String>[] idLists   = new HashSet[3];
        idLists[0] = lowGCIds; idLists[1] = midGCIds; idLists[2] = highGCIds;

        return idLists;
    }

    public static String[] correlationComparison(String comparisonGCT, String rnaSeqGCT,String OUT_DIR) throws IOException, InterruptedException {
        GCTFile rnaSeq = new GCTFile(rnaSeqGCT);
        String corrFiles[] = {OUT_DIR + "/corrRef.txt",OUT_DIR+"/corrMatrixSpearman.txt",OUT_DIR+"/corrMatrixPearson.txt"};

        if (comparisonGCT != null){
            //load comparison CGT file:
            GCTFile gct = new GCTFile(comparisonGCT);

            GCTFile combined  = rnaSeq.combine(gct, false, false); // combine by ID and not description

            String expFile =OUT_DIR + "/expressionForCorr.gct";
            combined.toFile(expFile);  // output combined file

            // run R script
            String referenceSamp = gct.getSamples()[0];
            runCorrelation(referenceSamp, combined, corrFiles[0]);
        }

        if (rnaSeq.getSamples().length > 1) {
            runCorrelationMatrix(rnaSeq, corrFiles[1]);
            runCorrelationMatrix(rnaSeq, corrFiles[2]);
        }
        
        return corrFiles;
    }


    public static void runCorrelation(String referenceSample, GCTFile gct,  String outfile) throws IOException, InterruptedException {
        String[] samples = gct.getSamples();
        int refIndex = -1;
        int index =0;
        while(refIndex == -1 && index < samples.length) {
            if(samples[index].equals(referenceSample)) {
                refIndex = index;
            }
            index++;
        }

        if(refIndex == -1) {
            System.err.println("Could not find reference sample" + referenceSample);
            System.exit(1);
        }

        PrintWriter writer = null;
        try {
            // spearman
            double[]spearmans = spearman(gct,refIndex);
            double[] pearsons = pearson(gct,refIndex,false);
            
            writer = new PrintWriter(outfile);
            writer.write("Sample\tSpearman\tPearson\n");
            for (int i = 0; i < samples.length; i++){
                writer.write(samples[i]); writer.write('\t');
                writer.write(""+spearmans[i]); writer.write('\t');
                writer.write(""+pearsons[i]); writer.write('\n');
            }
        } catch(IOException io) {
            io.printStackTrace();
            System.err.println("Unable to perform expression correlation comparison.");
        } finally {
            if(writer != null) {
                writer.close();
            }
        }
    }


    public static void runCorrelationMatrix( GCTFile gct,  String outfile) throws IOException, InterruptedException {

        String[] samps = gct.getSamples();
        boolean isPearson = outfile.toLowerCase().contains("pearson");
        
        BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
        // prepare header
        for (String samp: samps){
            out.write("\t");
            out.write(samp);
        }
        out.write("\n");

        // make matrix
        for (int i=0 ; i < samps.length; i++ ) {
            double [] results;
            if (isPearson){
                results = pearson(gct,i,true);
            }else{
                results = spearman(gct, i);
            }
            out.write(samps[i]);// the sample for this row
            
            // iterate over all samples as columns:
            for (int j = 0; j< samps.length; j++){
                out.write('\t');
                out.write(""+results[j]);
            }
            out.write('\n');
        }
        out.close();
    }


    public static double[] spearman(GCTFile gct, int refIndex){

        SpearmansCorrelation corr = new SpearmansCorrelation();
        double[] refColumn = gct.getColumnData(refIndex);

        ArrayList<Integer> indices = new ArrayList();
        for(int p =0;p <refColumn.length;p++)
        {
            if(!Double.isNaN(refColumn[p]))
            {
                indices.add(p);
            }
        }

        double[] refColTemp = new double[indices.size()];
        int nSamps = gct.getSamples().length;
        double[] results = new double[nSamps];
        for(int c=0; c < nSamps; c++) {
            double[] column = gct.getColumnData(c);

            double[] columnTemp = new double[indices.size()];

            //match NAs in both data columns
            for(int q =0;q < indices.size();q++) {
                columnTemp[q] = column[indices.get(q)];
                refColTemp[q] = refColumn[indices.get(q)];
            }

            double result = corr.correlation(refColTemp, columnTemp);
            results[c] = result;
        }
        return results;
    }


    /**
     * returns the Pearson Correlation between all samples in gct and the gct entry given by the refIndex.
     *
     * The data for each sample is logged, because it is RPKM, but the refData can be optionally logged depending
     * on whether it had been previously (e.g. normalized Affy)
     *
     * @param gct
     * @param refIndex
     * @param logRef
     * @return
     */
    public static double[] pearson(GCTFile gct, int refIndex, boolean logRef){

        PearsonsCorrelation corr = new PearsonsCorrelation();
        double[] refColumn = gct.getColumnData(refIndex);

        ArrayList<Integer> indices = new ArrayList();
        for(int p =0;p <refColumn.length;p++) {
            if(!Double.isNaN(refColumn[p]))
            {
                indices.add(p);
            }
        }

        double[] refColTemp = new double[indices.size()];
        int nSamps = gct.getSamples().length;
        double[] results = new double[nSamps];
        for(int c=0; c < nSamps; c++) {
            double[] column = gct.getColumnData(c);
            double[] columnTemp = new double[indices.size()];

            //match NAs in both data columns
            for(int q =0;q < indices.size();q++) {
                columnTemp[q] = Math.log(.01 + column[indices.get(q)]);
                refColTemp[q] = (logRef?  Math.log(0.01 +refColumn[indices.get(q)]):refColumn[indices.get(q)]);
            }
            double result = corr.correlation(refColTemp, columnTemp);
            results[c] = result;
        }
        return results;
    }


    public static int getTotReadsFromMetricFile(String metricFile) throws IOException{
        ArrayList<String[]> result= IOTools.fileToListofArrays(metricFile, "\\t");
        if (result == null) throw new FileNotFoundException("Couldn't find rRNA metrics file: " + metricFile);
        return Integer.parseInt(result.get(1)[0]);
    }


    /**
     * Turns the correlation file into html
     *
     * @param samples
     * @param corFile
     * @return
     * @throws IOException
     */
    public String createCorrHTMLTable(ArrayList<MetricSample> samples,String corFile) throws IOException {
        StringBuilder table = new StringBuilder();
        ArrayList<String[]> rows = IOTools.fileToListofArrays(corFile, null);
        ArrayList<String[]> vals = new ArrayList<String[]>();
        for (int i = 0 ; i < samples.size(); i++){ // using the samples.length because the vals.size will have an extra entry
            String[] vArray = new String[2];
            // skipping the first column, which is a repeat of the sample identifier
            vArray[0]=rows.get(i+1)[1]; // plus 1 to skip the header row
            vArray[1]=rows.get(i+1)[2];
            vals.add(vArray);
        }
        String[] header = new String[2]; header[0] ="Spearman"; header[1] ="Pearson";

        table.append("\n<h3>Correlation to Reference Expression Profile</h3>\n");
        table.append( getMetricTable(samples, vals, header,"11", metricTracker));
        table.append("\n<p><b>Spearman</b> is a ranked correlation coefficient.<b>Pearson</b> is correlation based on value.</p>\n");
        return table.toString();
    }


    /**
     * Turns the correlation file into html
     *
     * @param samples
     * @param corFile
     * @return
     * @throws IOException
     */
    public String createCorrHTMLTableMatrix(ArrayList<MetricSample> samples,String corFile) throws IOException {
        StringBuilder table = new StringBuilder();
        ArrayList<String[]> rows = IOTools.fileToListofArrays(corFile, null);
        String[] headerRow = rows.get(0);
        ArrayList<String[]> vals = new ArrayList<String[]>();
        for (int i = 0 ; i < samples.size(); i++){ // using the samples.length because the vals.size will have an extra entry
            String[] vArray = new String[headerRow.length-1];
            for (int j = 0 ; j < headerRow.length-1; j++){
                vArray[j]=rows.get(i+1)[j+1]; // i + 1 to skip the first row, which is the header; j+1 to skip the first col, which is sample id
            }            
            vals.add(vArray);
        }
        String[] header = new String[headerRow.length-1];
        String digitString = "";
        for (int i = 1 ; i < headerRow.length; i++){
            header[i-1] = headerRow[i];
            digitString+="1";
        }
        table.append( getMetricTable(samples, vals, header,digitString, metricTracker));
        return table.toString();
    }


    private static ArrayList<String> getOrderedContigsNoIndexing(String ref_genome) throws IOException {
        ArrayList<String> contigs = new ArrayList<String>();

        BufferedReader in = new BufferedReader(new FileReader(ref_genome)); // todo change this to work from the fast.fai index file
        String line = in.readLine();
        while(line !=null){
            try{
                if (line.startsWith(">")){
                    StringTokenizer toks = new StringTokenizer((line.substring(1)));
                    String contig = toks.nextToken();
                    contigs.add(contig);
                }
                line = in.readLine();
            }catch (RuntimeException e){
                System.out.println("Error parsing lines looking for contig names in reference");
                System.out.println("Line:\t" + line);
                e.printStackTrace();
                throw e;
            }
        }
        in.close();
        return contigs;
    }


    private static ArrayList<String> getOrderedContigsFAI(String ref_genome) throws IOException {
        ArrayList<String> contigs = new ArrayList<String>();

        BufferedReader in = new BufferedReader(new FileReader(ref_genome+".fai")); // todo change this to work from the fast.fai index file
        String line = in.readLine();
        while(line !=null){
            try{
                String contig = line.substring(0,line.indexOf('\t')); // contig name is hte first field
                contigs.add(contig);

                line = in.readLine();
            }catch (RuntimeException e){
                System.out.println("Error parsing lines looking for contig names in reference fai file: " + ref_genome+".fai");
                System.out.println("Line:\t" + line);
                e.printStackTrace();
                throw e;
            }
        }
        in.close();
        return contigs;
    }


    public static String getMetricTable(ArrayList<MetricSample> sampleInfo,
                                        ArrayList<String[]> metricResults, String[] header,
                                        String isDecimal, ArrayList<HashMap<String, String>> metricTracker) {
        try {
            StringBuilder str = new StringBuilder();

            str.append("<table><tr><th>Sample</th><th>Note</th>\n");
            for (String h: header){
                str.append("<th>").append(h).append("</th>");
            }
            str.append("\n</tr>\n");

            for (int i = 0 ; i < metricResults.size();i++){
                str.append("<tr>");
                str.append("<td>").append(sampleInfo.get(i).sampId).append("</td><td>").append(sampleInfo.get(i).notes).append("</td>");
                if (metricResults.get(i).length > isDecimal.length()){
                    System.out.println("Incompatible data and formatting");
                    System.out.println("\tHeader line " + Arrays.toString(header));
                    System.out.println("\tResults line " + Arrays.toString(metricResults.get(i)));
                    System.out.println("\tFormat line " + isDecimal);
                }
                for (int j = 0 ; j < metricResults.get(i).length; j++){

                    String cell = metricResults.get(i)[j];
                    metricTracker.get(i).put(header[j],cell);
                    String formatedNumber;
                    if (cell.equals("NA")){
                        formatedNumber = cell;
                    }else{
                        if (isDecimal.charAt(j) == '1'){
                            float f = Float.valueOf(cell);
                            if (Float.isNaN(f)|| f == 0f){
                                formatedNumber = "NA";
                            }else{
                                formatedNumber = String.format("%.3f", f);
                            }
                        }else{
                            int integer = Integer.valueOf(cell);
                            if (integer == 0){
                                formatedNumber = "NA";
                            }else{
                                formatedNumber = String.format("%,d", integer);
                            }
                        }
                    }
                    str.append("<td align='right'>").append(formatedNumber).append("</td>");
                }
                str.append("</tr>\n");
            }
            str.append("</table>\n");

            return str.toString();
        }catch(NumberFormatException e){
            System.err.println("Wrong number format for HTML Table");
            System.err.println("Format String: " + isDecimal);
            System.err.println("Header: " + Arrays.toString(header));
            throw e;
        }
    }

    public static String getStyle(){
        return "<style type=\"text/css\">\n" +
                "table {\n" +
                "\tborder-width: 1px;\n" +
                "\tborder-spacing: 2px;\n" +
                "\tborder-style: outset;\n" +
                "\tborder-color: gray;\n" +
                "\tborder-collapse: collapse;\n" +
                "\tbackground-color: rgb(209, 255, 255);\n" +
                "}\n" +
                "table th {\n" +
                "\tborder-width: 2px;\n" +
                "\tpadding: 2px;\n" +
                "\tborder-style: solid;\n" +
                "\tborder-color: gray;\n" +
                "\tbackground-color: rgb(209, 255, 255);\n" +
                "}\n" +
                "table td {\n" +
                "\tborder-width: 2px;\n" +
                "\tpadding: 2px;\n" +
                "\tborder-style: solid;\n" +
                "\tborder-color: gray;\n" +
                "\tbackground-color: white;\n" +
                "}\n" +
                "</style>"; //todo add style for <p> tag 'tLegend' (table legend)
    }


    private void createTopLevelIndex(ArrayList<MetricSample> samps, int MAX_TRANSCIRPTS,
                                            String metricTable, String corTable, String endLength,String outfile) throws IOException{
        BufferedWriter out = new BufferedWriter (new FileWriter(outfile));

        out.write("<html>\n<head><title>RNA-seq Metrics</title>\n");
        out.write(getStyle());
        out.write("\n</head>\n<body>\n");

        if (metricTable != null){
            out.write("<h1>RNA-seq Metrics</h1>\n");
            out.write("<h2>Read Count Metrics</h1>\n");

            out.write("<p>\n The following summary statistics are calculated by counting the number of reads that have the given characteristics.\n</p>\n");
            out.write(metricTable);
            out.write("\n");
        }

        if (corTable !=null) out.write(corTable); // correlation table

        out.write("\n<HR WIDTH=\"100%\" COLOR=\"#9090900\" SIZE=\"3\">\n");
        outputDoCResult(samps, out,"Coverage Metrics for Bottom "+MAX_TRANSCIRPTS+" Expressed Transcripts",PerBaseDoC.ExpressionLevel.low,"lowexpr", endLength, "bottom " + MAX_TRANSCIRPTS, null);
        outputDoCResult(samps, out,"Coverage Metrics for Middle "+MAX_TRANSCIRPTS+" Expressed Transcripts",PerBaseDoC.ExpressionLevel.medium,"medexpr", endLength, "middle " + MAX_TRANSCIRPTS, metricTracker);
        outputDoCResult(samps, out,"Coverage Metrics for Top "+MAX_TRANSCIRPTS+" Expressed Transcripts",PerBaseDoC.ExpressionLevel.high,"highexpr",  endLength, "top " + MAX_TRANSCIRPTS, null);
        out.write("\n<HR WIDTH=\"100%\" COLOR=\"#9090900\" SIZE=\"3\">\n");

        out.write("<h3>Mean Coverage</h3>\n");
        out.write("<p>" +
                "<h4>Low Expressed</h4>\n" +
                "<img src= 'meanCoverageNorm_low.png'>" +

                "<h4>Medium Expressed</h4>\n" +
                "<img src= 'meanCoverageNorm_medium.png'>" +

                "<h4>High Expressed</h4>\n" +
                "<img src= 'meanCoverageNorm_high.png'>"+
                "</p>\n");

        out.write("<h3>Mean Coverage from 3' End</h3>\n");
        out.write("<p>" +
                "<h4>Low Expressed</h4>\n" +
                "<img src= 'meanCoverage_low.png'>" +

                "<h4>Medium Expressed</h4>\n" +
                "<img src= 'meanCoverage_medium.png'>" +

                "<h4>High Expressed</h4>\n" +
                "<img src= 'meanCoverage_high.png'>"+
                "</p>\n");

        out.write("\n<HR WIDTH=\"100%\" COLOR=\"#9090900\" SIZE=\"3\">\n");
        
        if (gapLengthDistribution){
            out.write("<h3>Gap Length Distributions</h3>\n");
            out.write("<p>" +
                "<h4>Low Expressed</h4>\n" +
                "<img src= 'gapLengthHist_low.png'>" +

                "<h4>Medium Expressed</h4>\n" +
                "<img src= 'gapLengthHist_medium.png'>" +

                "<h4>High Expressed</h4>\n" +
                "<img src= 'gapLengthHist_high.png'>"+
                "</p>\n");
        }

        if (this.doStrat) {  // a stratification by GC analysis was performed
            out.write("\n<HR WIDTH=\"100%\" COLOR=\"#9090900\" SIZE=\"3\">\n");
            out.write("<p>\n");
            out.write("<h2>GC Stratification</h2>\n");
            out.write("<ul>\n");
            out.write("\t<li><a target='_new' href='gc/high/index.html'>High GC</a></li>\n");
            out.write("\t<li><a target='_new' href='gc/mid/index.html'>Moderate GC</a></li>\n");
            out.write("\t<li><a target='_new' href='gc/low/index.html'>Low GC</a></li>\n");
            out.write("</ul>");
            out.write("</p>\n");
        }

        out.write("\n<HR WIDTH=\"100%\" COLOR=\"#9090900\" SIZE=\"3\">\n");
        out.write("<p>\n");
        out.write("<h2>Files</h2>\n");
        out.write("<table><tr><th>File</th><th>Description</th></tr>\n");
        File curDir = new File("./");
        File outputfileDir = new File(outfile);
        outputfileDir = outputfileDir.getParentFile();
        out.write("<tr><td><a target='_new' href='metrics.tsv'>Metrics Tab Separated Value File</a></td><td>Text file containing all the metrics of the report in a single tab delimited file.</td></tr>\n");

        if (curDir.getName().equals(outputfileDir.getName())) {
            out.write("<tr><td><a target='_new' href='genes.rpkm.gct'>RPKM Values</a></td><td>A GCT file containing the expression profiles of each sample</td></tr>\n");
            out.write("<tr><td><a target='_new' href='countMetrics.html'>Read Count Metrics</a></td><td>An HTML file containing only the read count-based metrics</td></tr>\n");
        }

        out.write("<tr><td><a target='_new' href='meanCoverageNorm_low.txt'>Mean Coverage Plot Data - Low Expr</a></td><td>Text file containing the data for mean coverage plot by position for low expression coverage</td></tr>\n");
        out.write("<tr><td><a target='_new' href='meanCoverageNorm_medium.txt'>Mean Coverage Plot Data - Medium Expr</a></td><td>Text file containing the data for mean coverage plot by position for medium expression coverage</td></tr>\n");
        out.write("<tr><td><a target='_new' href='meanCoverageNorm_high.txt'>Mean Coverage Plot Data - High Expr</a></td><td>Text file containing the data for mean coverage plot by position for high expression coverage</td></tr>\n");

        out.write("<tr><td><a target='_new' href='meanCoverage_low.txt'>Mean Coverage Plot Data - Low Expr</a></td><td>Text file containing the data for mean coverage plot by distance from 3' end for low expression coverage</td></tr>\n");
        out.write("<tr><td><a target='_new' href='meanCoverage_medium.txt'>Mean Coverage Plot Data - Medium Expr</a></td><td>Text file containing the data for mean coverage plot by distance from 3' end for medium expression coverage</td></tr>\n");
        out.write("<tr><td><a target='_new' href='meanCoverage_high.txt'>Mean Coverage Plot Data - High Expr</a></td><td>Text file containing the data for mean coverage plot by distance from 3' end for high expression coverage</td></tr>\n");

        out.write("<tr><td><a target='_new' href='gapLengthHist_low.txt'>Mean Coverage Plot Data - Low Expr</a></td><td>Text file containing the data for gap length counts for low expression coverage</td></tr>\n");
        out.write("<tr><td><a target='_new' href='gapLengthHist_medium.txt'>Mean Coverage Plot Data - Medium Expr</a></td><td>Text file containing the data for gap length counts for medium expression coverage</td></tr>\n");
        out.write("<tr><td><a target='_new' href='gapLengthHist_high.txt'>Mean Coverage Plot Data - High Expr</a></td><td>Text file containing the data for gap length counts for high expression coverage</td></tr>\n");
        out.write("\n</table>\n");

        out.write("</p>\n");

        out.write("\n<HR WIDTH=\"100%\" COLOR=\"#9090900\" SIZE=\"3\">\n");
        out.write("<p>\n");
        out.write("<h2>Summary of Runtime Parameters</h2>\n");
        out.write("<table><tr><th>Option</th><th>Description</th><th>Value</th></tr>\n");
        out.write("<tr><td>Samples</td><td>Samples/Sample File used</td><td>");
        out.write(Tools.stripFileName(this.SAMPLE_FILE)); out.write("</td></tr>\n");
        out.write("<tr><td>Transcript Model</td><td>GTF formatted file containing the transcript definitions</td><td>");
        out.write(Tools.stripFileName(this.TRANSCRIPT_MODEL)); out.write("</td></tr>\n");
        out.write("<tr><td>Reference Genome</td><td>The genome version to which the BAM is aligned</td><td>");
        out.write(Tools.stripFileName(this.REF_GENOME)); out.write("</td></tr>\n");
        out.write("<tr><td>Downsampling</td><td>For Coverage Metrics, the number of reads is randomly reduced to the given level</td><td>");
        out.write(this.dSampleTarget!=null?this.dSampleTarget:"none"); out.write("</td></tr>\n");
        out.write("<tr><td>Detailed Report</td><td>The optional detailed report contains coverage metrics for every transcript</td><td>");
        out.write(this.details?"details included":"no details"); out.write("</td></tr>\n");
        out.write("<tr><td>rRNA Intervals</td><td>Genomic coordinates of rRNA loci</td><td>");
        out.write(this.rRNAFile!=null?Tools.stripFileName(this.rRNAFile):"taken from GTF file"); out.write("</td></tr>\n");
        out.write("\n</table>\n");

        out.write("</p>\n");
        out.write("<h5>Generated by <a target='_newWebPage' href='http://www.broadinstitute.org/rna-seqc'>RNA-SeQC</a> "+VERSION+"</h5>\n");
        out.write("<h5>Run on "+new Date()+"</h5>\n");
        out.write("</body>\n</html>");

        out.close();
    }

    private static void outputDoCResult(ArrayList<MetricSample> samples, BufferedWriter out,
                                        String title, PerBaseDoC.ExpressionLevel level, String dir, String endLength, String footnote,
                                        ArrayList<HashMap<String, String>> metricTracker) throws IOException {
        out.write("<h2>"+title+"</h1>\n");
        out.write("<p>\n The metrics in this table are calculated across the transcripts that were determined to have the highest expression levels. ");

        out.write("<table border=1>");
        out.write("<tr><th>Sample</th><th>Note</th><th>Mean Per Base Cov.</th><th>Mean CV</th>");
        out.write("<th>No. Covered 5'</th><th>5'"+endLength+"Base Norm</th>");
        out.write("<th>No. Covered 3'</th><th>3' "+endLength+"Base Norm</th>");
        out.write("<th>Num. Gaps</th><th>Cumul. Gap Length</th><th>Gap %</th>");

        out.write("</tr>\n");
        int i = 0;
        for (MetricSample samp: samples){
            out.write("<tr>");
            out.write("<td><a href='"+samp.sampId+"/"+dir +"/index.html'>"+samp.sampId + "</a></td>");
            out.write("<td>"+samp.notes + "</td>");

            //vvout.write("<td align='right'>"+doc.getTotalCoverage()+"</td>");
            //todo ADD STANDARD DEVIATIONS TO ALL THESE CALCULATIONS

            TranscriptListDoCResult doc = new TranscriptListDoCResult();
            doc.loadSummaryFromFile(samp.getDoCTranscriptsSummaryResultsFile(level));

            out.write("<td align='right'>"+String.format("%.2f",doc.getAverageCoverages())+"</td>");
            if (metricTracker !=null) metricTracker.get(i).put("Mean Per Base Cov.",""+doc.getAverageCoverages());
            out.write("<td align='right'>" + String.format("%.2f", doc.getAverageCV()) + "</td>");
            if (metricTracker !=null) metricTracker.get(i).put("Mean CV",""+doc.getAverageCV());
            out.write("<td align='right'>"+ doc.getFiveEndCovCount() + "</td>"); //toodo add endLength as param to this as well
            if (metricTracker !=null) metricTracker.get(i).put("No. Covered 5'",""+doc.getFiveEndCovCount());

            if (endLength.equals("200")){
                out.write("<td align='right'>"+ String.format("%.2f",doc.getAvFiveEnd200()) + "</td>");
                if (metricTracker !=null) metricTracker.get(i).put("5' Norm",""+doc.getAvFiveEnd200());
            }else if(endLength.equals("50")){
                out.write("<td align='right'>"+String.format("%.2f", doc.getAvFiveEnd50())+ "</td>");
                if (metricTracker !=null) metricTracker.get(i).put("5' Norm",""+doc.getAvFiveEnd50());
            }else if (endLength.equals("100")){
                out.write("<td align='right'>"+String.format("%.2f", doc.getAvFiveEnd100()) + "</td>");
                if (metricTracker !=null) metricTracker.get(i).put("5' Norm",""+doc.getAvFiveEnd100());
            }

            out.write("<td align='right'>"+ doc.getThreeEndCovCount() + "</td>");
            if (endLength.equals("200")){
                out.write("<td align='right'>"+ String.format("%.3f",doc.getAvThreeEnd200()) + "</td>");
                if (metricTracker !=null) metricTracker.get(i).put("3' Norm",""+doc.getAvThreeEnd200());
            }else if(endLength.equals("50")){
                out.write("<td align='right'>"+String.format("%.3f",doc.getAvThreeEnd50()) + "</td>");
                if (metricTracker !=null) metricTracker.get(i).put("3' Norm",""+doc.getAvThreeEnd50());
            }else if (endLength.equals("100")){
                out.write("<td align='right'>"+String.format("%.3f", doc.getAvThreeEnd100()) + "</td>");
                if (metricTracker !=null) metricTracker.get(i).put("3' Norm",""+doc.getAvThreeEnd100());
            }

            out.write("<td align='right'>"+ doc.getNumberOfGaps() + "</td>");
            if (metricTracker !=null) metricTracker.get(i).put("Num. Gaps",""+doc.getNumberOfGaps());
            out.write("<td align='right'>"+doc.getCumulativeGapLength() + "</td>");
            if (metricTracker !=null) metricTracker.get(i).put("Cumul. Gap Length",""+doc.getCumulativeGapLength());
            out.write("<td align='right'>"+String.format("%.1f",(doc.getGapRatio() * 100f)) + "</td>");
            if (metricTracker !=null) metricTracker.get(i).put("Gap %",""+doc.getGapRatio());

            out.write("</tr>\n");
            i++;
        }
        out.write("</table>\n");

        out.write("<p>\n");
        out.write("It is important to note that these values are restricted to the "+footnote+" expressed transcripts. ");
        out.write(" <b>5'</b> and <b>3'</b> values are per-base coverage averaged across all top transcripts.  ");
        out.write(" 5' and 3' ends are "+endLength+" base pairs. ");
        out.write(" Gap % is the total cumulative gap length divided by the total cumulative transcript lengths.");
        out.write("\n</p>\n");
    }

    private static String getSummaryRow(String summaryFile) {
        try{
            BufferedReader in = new BufferedReader (new FileReader(summaryFile.trim()));
            in.readLine(); // skip header;
            String result = "";
            String line = in.readLine(); // first and only line

            for (String cell: line.split("\\t")) {
                result+="<td align='right'>"+cell+"</td>";
            }

            in.close();
            return result;

        } catch(IOException e) {
            e.printStackTrace();
            return "<tr><td>could not find summary file</td></tr>";
        }
    }


    /**
     *
     * @param TRANSCRIPT_GTF
     * @param REF_GENOME
     * @param refGeneFile
     * @return
     * @throws IOException
     */
    public static String createRefGeneAndRRNAFiles(String TRANSCRIPT_GTF, String REF_GENOME, String refGeneFile, String rRNAIntervalFile, String transcriptTypeField) throws IOException {
        System.out.println("Retriving contig names from reference");
        ArrayList<String> contigs = getOrderedContigsFAI(REF_GENOME);
        System.out.println("\t contig names in reference: "+contigs.size());


        System.out.println("Loading GTF for Read Counting");
        TranscriptList transcripts = TranscriptList.loadGTF(TRANSCRIPT_GTF, transcriptTypeField);
        System.out.println("Converting to refGene");

        Transcript.transcriptsToRefGeneFile(contigs, transcripts, refGeneFile);
        if (rRNAIntervalFile != null){
            transcripts.toRRNAIntervalList(rRNAIntervalFile);
        }

        return refGeneFile;
    }


    public static class MetricSample {
        String sampId;
        String bamFile;
        String notes;
        String listFile = null;
        private String tmpMetricFile;
        private String sampDir;
        private String rRNAMetricsFile;
        private String libraryComplexityFile;
        private String metricsFile;
        private String gctFile;

        public MetricSample(String samp, String bam, String  note) {
            sampId = samp;
            bamFile = bam;
            notes = note;
            listFile = null;
        }

        /**
         * Reads in the sample names, bam files paths and notes for each sample
         *
         * To support multiple bams per sample, the parser will look for multiple instances of the sample ID
         * If more than one instance is found, a new List file for GATK is created, the sample is condensed, and the
         * list file replaces the bam file
         *
         *
         * @param sampleFile
         * @return
         * @throws IOException
         */
        public static ArrayList<MetricSample> readInSamples(String sampleFile, String outDir) throws IOException{
            if (sampleFile.contains("|")) {
                return readInSamplesFromCL(sampleFile,outDir,"\\|");
            } else if (sampleFile.contains(",")) {
                return readInSamplesFromCL(sampleFile,outDir,",");
            }

            HashMap<String,ArrayList<String[]>> samples = new HashMap<String,ArrayList<String[]>>();

            BufferedReader in = new BufferedReader(new FileReader(sampleFile));
            String line =in.readLine();
            line = in.readLine(); // skip header
            ArrayList<String> origSampOrder = new ArrayList<String>();
            while (line!=null) {
                if (!line.trim().startsWith("#") && !line.trim().equals("")) {
                    String[] split = line.split("\\t");
                    //trim
                    for (int i = 0 ; i < split.length; i++) {
                        if (split[i]!=null){
                            split[i] = split[i].trim();
                        }
                    }
                    ArrayList<String[]> sampList = samples.get(split[0]);
                    if (sampList == null) {
                        sampList = new ArrayList<String[]>();
                        samples.put(split[0],sampList);
                        origSampOrder.add(split[0]);
                    }
                    sampList.add(split);
                }
                line = in.readLine();
            }

            ArrayList<MetricSample> finalSamps = new ArrayList<MetricSample>(origSampOrder.size());
            for (String s: origSampOrder) {
                ArrayList<String[]> thisList = samples.get(s);
                if (thisList.size() == 1) {
                    String[] samp = thisList.get(0);
                    finalSamps.add(new MetricSample(samp[0],samp[1],samp[2]));
                } else {
                    // make list file
                    String bamList = "";
                    for (String[] split: thisList) {
                        bamList+=split[1]+"\n";
                    }
                    String listFile = outDir+ "/" + s+".list";
                    IOTools.stringToFile(listFile,bamList);
                    String[] l = thisList.get(0);
                    MetricSample ms = new MetricSample(l[0],l[1],l[2]);
                    ms.setListFile(listFile);
                    finalSamps.add(ms);
                }
            }
            MetricSample.prepare(finalSamps, outDir);
            return finalSamps;
        }


        private static ArrayList<MetricSample> readInSamplesFromCL(String sampleString, String outDir, String delim) {
            ArrayList<MetricSample> samps = new ArrayList<MetricSample>();
            String[] l = sampleString.split(delim);
            MetricSample ms = new MetricSample(l[0],l[1],l[2]);
            samps.add(ms);
            MetricSample.prepare(samps, outDir);
            return samps;
        }

        /**
         * creates the sample directories for all samples
         * names the temporary sample file for all samples
         */
        private static void prepare(ArrayList<MetricSample> samples, String outDir) {
            for (RNASeqMetrics.MetricSample samp: samples) {
                String sampDir = outDir + "/" + samp.sampId;
                File dir = new File(sampDir);
                dir.mkdir();

                String tmpMetricFile =   sampDir +"/" + samp.sampId + ".metrics.tmp.txt";
                String rRNAMetricsFile = sampDir +"/" + samp.sampId + ".rRNA_counts.txt";
                String lcOutputFile =    sampDir +"/" + samp.sampId + ".libraryComplexity.txt";
                String metricsFile =     sampDir +"/" + samp.sampId + ".metrics.txt";
                String gctFile =         sampDir +"/" + samp.sampId + ".transcripts.rpkm.gct";

                samp.setSampDir(sampDir);
                samp.setTmpMetricFile(tmpMetricFile);
                samp.setrRNAMetricsFile(rRNAMetricsFile);
                samp.setLibraryComplexityFile(lcOutputFile);
                samp.setMetricsFile(metricsFile);
                samp.setGCTFile(gctFile);
            }
        }

        public void setListFile(String listFile) {
            this.listFile = listFile;
        }

        public boolean hasList() {
            return this.listFile != null;
        }

        public void setTmpMetricFile(String tmpMetricFile) {
            this.tmpMetricFile = tmpMetricFile;
        }

        public void setSampDir(String sampDir) {
            this.sampDir = sampDir;
        }

        public String getBamFileOrList() {
            if (this.hasList()) {
                return this.listFile;
            } else {
                return this.bamFile;
            }
        }

        public String getTmpMetricsFile() {
            return tmpMetricFile;
        }

        public void setrRNAMetricsFile(String rRNAMetricsFile) {
            this.rRNAMetricsFile = rRNAMetricsFile;
        }

        public String getrRNAMetricsFile() {
            return rRNAMetricsFile;
        }

        public String getSampleDirectory() {
            return this.sampDir;
        }

        public void setLibraryComplexityFile(String libraryComplexityFile) {
            this.libraryComplexityFile = libraryComplexityFile;
        }

        public String getLibraryComplexityFile() {
            return libraryComplexityFile;
        }

        public void setMetricsFile(String metricsFile) {
            this.metricsFile = metricsFile;
        }

        public String getMetricsFile() {
            return this.metricsFile;
        }

        public String getGCTFile() {
            return gctFile;
        }

        public void setGCTFile(String file){
            this.gctFile = file;
        }

        public int getTotalReads() throws IOException {
            String metricFile = this.getTmpMetricsFile();
            int totReads = getTotReadsFromMetricFile(metricFile);
            return totReads;
        }

        public String getExpressionList(PerBaseDoC.ExpressionLevel lev) {
            return getExpressionDir(lev)+sampId+ PerBaseDoC.EXPR_SUFFIX;
        }

        public String getExpressionDir(PerBaseDoC.ExpressionLevel level) {
            switch(level){
                case low: return sampDir+"/lowexpr/";
                case medium: return sampDir+"/medexpr/";
                case high: return sampDir+"/highexpr/";
                // case all: return sampDir+"/allexpr/";  // breaks pre-1.1.9 version
            }
            return null;
        }

        public String getIntervalFile(PerBaseDoC.ExpressionLevel level) {
            return getExpressionDir(level)+"intervals.list";
        }

        public String getDoCResultsFile(PerBaseDoC.ExpressionLevel level) {
            return getExpressionDir(level) +"/perBaseDoC.out";
        }

        public String getMeanCoverageByPosFile(PerBaseDoC.ExpressionLevel level) {
            return getExpressionDir(level)+"/meanCovByPosition.txt";
        }

        public String getMeanCoverageByPosNormalizedFile(PerBaseDoC.ExpressionLevel level) {
            return getExpressionDir(level)+"/meanCovByPositionNormed.txt";
        }

        public String getGapLengthHistogramFile(PerBaseDoC.ExpressionLevel level) {
            return getExpressionDir(level)+"/gapLengthHistogram.txt";
        }

        public String getDoCTranscriptsResultsFile(PerBaseDoC.ExpressionLevel level) {
            return getExpressionDir(level)+sampId+".DoCTranscripts";
        }

        public String getDoCTranscriptsSummaryResultsFile(PerBaseDoC.ExpressionLevel level) {
            return getExpressionDir(level)+sampId+".DoCTranscriptsSummary";
        }
    }

}
