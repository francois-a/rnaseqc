package org.broadinstitute.cga.rnaseq.gatk;

import org.broadinstitute.cga.tools.Performance;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: ddeluca
 * Date: 3/31/11
 * Time: 10:18 AM
 * To change this template use File | Settings | File Templates.
 */
public class GATKTools {

    public static final boolean QUIET = true; //todo quiet is gone from GATK. where did it go

    public static void main(String[] args) {
        try{
            System.out.println("REF_GENOME BAM_FILE INTERVAL REFSEQ  outfile");
            runIntronReadCountDebug(args[0],args[1],args[2],args[3],args[4]);
        }catch(Exception e){
            e.printStackTrace();
        }
    }

    /**
     * Some other parameters to consider supporting:
     *
     *  -stddev
     *  -mmq 30 -mbq 20
     *
     *
     * @param REF_GENOME
     * @param BAM_FILE
     * @param docResults
     * @param intervalList
     * @param ommitBases
     * @param extraFlags example: [no space]-rf Downsampling -numReads 2000 -targetReads 1000
     * @throws Exception
     */
    public static void runDoC(String REF_GENOME, String BAM_FILE, String docResults,String intervalList, boolean ommitBases, String extraFlags) throws Exception {
        // example:java -jar ../../scripts/GenomeAnalysisTK.jar -T DepthOfCoverage -R /xchip/cga2/berger/RNASeq/refs/Ensembl52.plus.Genome.fasta -I ../CLL-JP/CLL-JP.rna.GATKRecalibrated.flagged.bam -o ./myDoC.out -L top_expressed_intervals.list
        
        System.out.println("Running GATK Depth of Coverage Analysis ....");
        Performance perf = new Performance("Depth of Coverage run time",Performance.Resolution.minutes);
        CommandLineGATK instance = new CommandLineGATK();
        String commandStr = "-T DepthOfCoverage"+  (extraFlags!=null?" "+extraFlags:"")+" -R "+REF_GENOME+" -I "+BAM_FILE+" -o "+docResults+
                " -L "+intervalList + (ommitBases?" -omitBaseOutput":"") +(QUIET?" -l ERROR":"") ;

        String argv[] = commandStr.split(" ");
        System.out.println("Arguments:\t"+ commandStr);
        System.out.println("Arguments Array:\t"+ Arrays.toString(argv));
        CommandLineGATK.start(instance, argv);

        System.out.println("GATK command result code: " + CommandLineProgram.result);
        System.out.println(perf);
        System.out.println("\t ... GATK Depth of Coverage Analysis DONE");
    }


    /**
     * Runs the CountReadMetrics walker
     * 
     * @param REF_GENOME
     * @param BAM_FILE
     * @param intervalList
     * @param refseq
     * @param outfile
     * @throws Exception
     */
    public static void runCountReadMetrics(String REF_GENOME, String BAM_FILE, String intervalList, String refseq, String outfile/*, boolean dSample*/) throws Exception {
        // example:java -jar ../../scripts/GenomeAnalysisTK.jar -T DepthOfCoverage -R /xchip/cga2/berger/RNASeq/refs/Ensembl52.plus.Genome.fasta -I ../CLL-JP/CLL-JP.rna.GATKRecalibrated.flagged.bam -o ./myDoC.out -L top_expressed_intervals.list
        System.out.println("Running CoutReadMetrics Walker ....");
        CommandLineGATK instance = new CommandLineGATK();

        String commandStr =  "-T CountReadMetrics --outfile_metrics "+outfile +
                /*(dSample?" -rf Downsampling -numReads 2000 -targetReads 1000":"" )+ */
                " -R "+REF_GENOME+" -I "+BAM_FILE+ (intervalList!=null?" -L "+intervalList:"") +
                //" -B:refseq,RefSeq "+refseq;
                " -refseq "+refseq+(QUIET?" -l ERROR":"");

        String argv[] = commandStr.split(" ");
        System.out.println("Arguments: "+ Arrays.toString(argv));
        CommandLineGATK.start(instance, argv);

        System.out.println("GATK command result code: " + CommandLineProgram.result);
        System.out.println("\t ... GATK CoutReadMetrics Analysis DONE");
    }

    /**
     * Runs the CountReadMetrics walker
     *
     * @param BAM_FILE
     * @param intervalList
     * @param refseq
     * @param outfile
     * @throws Exception
     */
    public static void runCountReadMetricsNoRef(String BAM_FILE, String intervalList, String refseq, String outfile/*, boolean dSample*/) throws Exception {
        System.out.println("Running IntronicExpressionReadBlock Walker ....");
              CommandLineGATK instance = new CommandLineGATK();

              String commandStr =  "-T IntronicExpressionReadBlock --outfile_metrics "+outfile +
                      //(dSample?" -rf Downsampling -numReads 2000 -targetReads 1000":"" ) +
                      " -I "+BAM_FILE+ (intervalList!=null?" -L "+intervalList:"") +
                      //" -B:refseq,RefSeq "+refseq;
                      " -refseq "+refseq+(QUIET?" -l ERROR":"");

              String argv[] = commandStr.split(" ");
              System.out.println("Arguments: "+ Arrays.toString(argv));
              CommandLineGATK.start(instance, argv);

              System.out.println("GATK command result code: " + CommandLineProgram.result);
              System.out.println("\t ... GATK CoutReadMetrics Analysis DONE");
    }

    /**
     * Runs the OutputCountReadsWalker walker
     * This simple counter just outputs the number of reads into a file
     *
     * @param REF_GENOME
     * @param BAM_FILE
     * @param intervalList
     * @param outfile
     * @throws Exception
     */
    public static void runIntervalReadCounter(String REF_GENOME, String BAM_FILE, String intervalList,
                                              String outfile, String additionalFlags/*, boolean dSample*/) throws Exception {
        // example:java -jar ../../scripts/GenomeAnalysisTK.jar -T DepthOfCoverage -R /xchip/cga2/berger/RNASeq/refs/Ensembl52.plus.Genome.fasta -I ../CLL-JP/CLL-JP.rna.GATKRecalibrated.flagged.bam -o ./myDoC.out -L top_expressed_intervals.list
        System.out.println("Running OutputCountReads Walker ....");
        CommandLineGATK instance = new CommandLineGATK();

        String commandStr =  "-T OutputCountReadsWalker --outfile_readCounts "+outfile +
                /*(dSample?" -rf Downsampling -numReads 2000 -targetReads 1000":"" )+ */
                " -R "+REF_GENOME+" -I "+BAM_FILE+ (intervalList!=null?" -L "+intervalList:"") +
                //" -B:refseq,RefSeq "+refseq;
                (QUIET?" -l ERROR":"")+ (additionalFlags!=null?" "+additionalFlags:"");

        String argv[] = commandStr.split(" ");
        System.out.println("Arguments: "+ Arrays.toString(argv));
        CommandLineGATK.start(instance, argv);

        System.out.println("GATK command result code: " + CommandLineProgram.result);
        System.out.println("\t ... GATK CoutReadMetrics Analysis DONE (file:"+outfile+")");
    }


    public static void test(String REF_GENOME, String BAM_FILE, String intervalList, String refseq, String outfile, boolean dSample) throws Exception {
        // example:java -jar ../../scripts/GenomeAnalysisTK.jar -T DepthOfCoverage -R /xchip/cga2/berger/RNASeq/refs/Ensembl52.plus.Genome.fasta -I ../CLL-JP/CLL-JP.rna.GATKRecalibrated.flagged.bam -o ./myDoC.out -L top_expressed_intervals.list
        System.out.println("Running TestReadWalker ....");
        CommandLineGATK instance = new CommandLineGATK();

        String commandStr =  "-T TestRead --outfile_metrics "+outfile +
                (dSample?" -rf Downsampling -numReads 2000 -targetReads 1000":"" )+
                " -R "+REF_GENOME+" -I "+BAM_FILE+ (intervalList!=null?" -L "+intervalList:"") +
                //" -B:refseq,RefSeq "+refseq;
                " -refseq "+refseq+(QUIET?" -l ERROR":"");

        String argv[] = commandStr.split(" ");
        System.out.println("Arguments: "+ Arrays.toString(argv));
        CommandLineGATK.start(instance, argv);

        System.out.println("GATK command result code: " + CommandLineProgram.result);
        System.out.println("\t ... GATK CoutReadMetrics Analysis DONE");
    }

    public static void test2(String REF_GENOME, String BAM_FILE) throws Exception {
        // example:java -jar ../../scripts/GenomeAnalysisTK.jar -T DepthOfCoverage -R /xchip/cga2/berger/RNASeq/refs/Ensembl52.plus.Genome.fasta -I ../CLL-JP/CLL-JP.rna.GATKRecalibrated.flagged.bam -o ./myDoC.out -L top_expressed_intervals.list
        System.out.println("Running TestReadWalker ....");

        System.out.println("This constructor should be empty..");
        CommandLineGATK instance = new CommandLineGATK();
        System.out.println("...constructed");

        String commandStr =  "-T TestRead"+
                " -R "+REF_GENOME+" -I "+BAM_FILE;

        String argv[] = commandStr.split(" ");
        System.out.println("Arguments: "+ Arrays.toString(argv));
        System.out.println("Starting");
        CommandLineGATK.start(instance, argv);

        System.out.println("GATK command result code: " + CommandLineProgram.result);
        System.out.println("\t ... GATK CoutReadMetrics Analysis DONE");
    }


    /**
     *
     * @param REF_GENOME
     * @param BAM_FILE
     * @param intervalList
     * @param refseq
     * @param outfile
     * @throws Exception
     */
    public static void runIntronReadCount(String REF_GENOME, String BAM_FILE, String intervalList, String refseq, String outfile, String additionalFlags,boolean strictMode) throws Exception {
        System.out.println("Running IntronicExpressionReadBlock Walker ....");
        CommandLineGATK instance = new CommandLineGATK();

        String commandStr =  "-T IntronicExpressionReadBlockWalker --outfile_metrics "+outfile +
                //(dSample?" -rf Downsampling -numReads 2000 -targetReads 1000":"" ) +
                " -R "+REF_GENOME+" -I "+BAM_FILE+ (intervalList!=null?" -L "+intervalList:"") +
                (strictMode?" -strict":"")+
                " -refseq "+refseq+(QUIET?" -l ERROR":"")+ (additionalFlags!=null?" "+additionalFlags:"");

        String argv[] = commandStr.split(" ");
        System.out.println("Arguments: "+ Arrays.toString(argv));
        CommandLineGATK.start(instance, argv);

        System.out.println("GATK command result code: " + CommandLineProgram.result);
        System.out.println("\t ... GATK CoutReadMetrics Analysis DONE");
    }

    /**
     *
     * @param REF_GENOME
     * @param BAM_FILE
     * @param intervalList
     * @param refseq
     * @param outfile
     * @throws Exception
     */
    public static void runIntronReadCountDebug(String REF_GENOME, String BAM_FILE, String intervalList, String refseq, String outfile) throws Exception {
        System.out.println("Running IntronicExpressionReadBlock Walker ....");
        CommandLineGATK instance = new CommandLineGATK();

        String commandStr =  "-T IntronicExpressionReadBlockDebug --outfile_metrics "+outfile +
                //(dSample?" -rf Downsampling -numReads 2000 -targetReads 1000":"" ) +
                " -R "+REF_GENOME+" -I "+BAM_FILE+ (intervalList!=null?" -L "+intervalList:"") +
                //" -B:refseq,RefSeq "+refseq;
                " -refseq "+refseq+(QUIET?" -l ERROR":"");

        String argv[] = commandStr.split(" ");
        System.out.println("Arguments: "+ Arrays.toString(argv));
        CommandLineGATK.start(instance, argv);

        System.out.println("GATK command result code: " + CommandLineProgram.result);
        System.out.println("\t ... GATK CoutReadMetrics Analysis DONE");
    }

}
