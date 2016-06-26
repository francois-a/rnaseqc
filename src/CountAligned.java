package org.broadinstitute.cga.picardbased;


import net.sf.picard.sam.MarkDuplicates;
import net.sf.picard.sam.SamToFastq;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.cga.tools.IOTools;
import org.broadinstitute.cga.tools.Tools;


import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * Created by the Cancer Genome Analysis Group at the Broad Institute.
 * Author: David S. DeLuca
 * User: ddeluca
 * Date: 7/19/11
 * Time: 10:55 AM
 * Developed as part of the RNA-seq analysis efforts of the Broad Institute
 */
public class CountAligned {
    public String BWA = "bwa" ; // default assumes bwa is in path


    private static final short DISTANCE_THRESHOLD = 5;
    private String reference;
    private String tmpDir ;
    private float[] result =null ;
    private int[] count;
    private boolean keepFiles = false;

    private double dSampleRate = 0d;


    public CountAligned(String reference, String tmpSpace) {
        this.reference = reference;
        this.tmpDir = tmpSpace;
    }

    public static void main(String[] args) {


        try{

            String reference = "/Users/ddeluca/Documents/RNA-seq/metrics/rRNA/human_all_rRNA.fasta";

            //"/Users/ddeluca/Documents/RNA-seq/metrics/Solexa-63857/B019RACXX.2.aligned.duplicates_marked.bam";//args[0];
            String tmpSpace =  "/Users/ddeluca/Documents/RNA-seq/metrics/rRNA/";

            String inputBam = tmpSpace + "nugen.bam";


            CountAligned ca = new CountAligned(reference,  tmpSpace);
            ca.BWA = "/Users/ddeluca/Documents/Code/Libraries/bwa-0.5.9/bwa";
            ca.keepFiles = true;

            ca.countBAM(inputBam, false);


            System.out.println("Result:" + ca.getRatio());

            System.out.println("Done.");
        }catch(Exception e){
            System.out.println("Fatal Error!");
            e.printStackTrace();
        }
    }


    /**
     * Performs:    1) samToFastq
     *              2) calls a fuction that runs the alignment and counts the reads that aligned
     * @param inputBam
     * @param singleEnd
     * @return
     * @throws IOException
     * @throws InterruptedException
     */
    public float[] countBAM(String inputBam, boolean singleEnd) throws IOException, InterruptedException {
        if (this.result == null){
            initResult();
        }


        //if reference is a sam file, then the rRNA has already been aligned
        if(reference!= null && reference.endsWith(".sam"))
        {
            countAligned(reference);
        }
        else
        {
            String out1 = tmpDir + "/e1.fastq";
            String out2 = null;

            if(!singleEnd)
            {
                out2 = tmpDir + "/e2.fastq";
            }


            String dSampledBam = null;
            if (dSampleRate > 0d && dSampleRate < 1d){
                System.out.println("Downsampling before aligning at rate: " + dSampleRate);
                dSampledBam = tmpDir +"/dSample.bam";
                String[]argv = {"I="+inputBam, "O="+dSampledBam, "P="+dSampleRate, "QUIET=true"};
                int rVal = (new org.broadinstitute.cga.picardbased.DownsampleSam()).instanceMain(argv);
                System.out.println("Downsampling exited with code: " + rVal);
                inputBam=dSampledBam;
            }


            System.out.println("Extracting FASTQ");
            getFastQ(inputBam,out1,out2);


            // this step runs the alignment as well as counts the reads that aligned
            countFastq(out1, out2);

            if(!keepFiles){
                new File(out1).delete();
                if (!singleEnd) new File(out2).delete();
                if(dSampledBam !=null)
                    new File(dSampledBam).delete();
            }
        }

        return result;

    }

    public float[] countFastq(String fastq1, String fastq2) throws IOException, InterruptedException {
        if (this.result == null){
            initResult();
        }

        System.out.println("BWA on end 1");
        String sai1 = tmpDir+"/end1.sai";
        runBWA( fastq1,sai1);

        String newSam = tmpDir+"/rRNA.sam";

        String sai2 = null;
        if(fastq2 != null)
        {
            System.out.println("BWA on end 2");
            sai2 = tmpDir+"/end2.sai";
            runBWA(fastq2, sai2);
            pairBWA( sai1,sai2,fastq1,fastq2,newSam);
        }
        else
        {
            singleBWA( sai1, fastq1,newSam);
        }


//        String samMarked =tmpDir+"/samMarked.sam";
//        String metrics  = tmpDir+"/picardMetrics.txt";
//        markDuplicates(newSam, samMarked, metrics);

        System.out.println("Counting aligned reads in " + newSam);

        countAligned(newSam);

        if(!keepFiles){
            new File(sai1).delete();
            if(fastq2 != null)
            {
                new File(sai2).delete();
            }
            new File(newSam).delete();
//            new File(samMarked).delete();
//            new File(metrics).delete();
        }

        return result;

    }

    private void markDuplicates(String inSam, String outSam, String metrics) {
        String[] args = {"I="+inSam, "F="+outSam, "METRICS_FILE="+metrics};

        MarkDuplicates.main(args);
    }

    private void initResult() {
        this.result = new float[6];
        Arrays.fill(result,0f);
    }

    private float[] countAligned(String inputFile) {

        //[0] = total reads ; [1] = aligned reads ; [2] = ratio; [3] = mapped with few mismatches [4] = ratio few mismatches; [5] aligned unique

        final SAMFileReader referenceSam = new SAMFileReader(new File(inputFile));

        referenceSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
//        int secs = 0;
//        int count = 0;
        int nmErrorLogCount = 0;
        for (final SAMRecord read : referenceSam) {

            if (!read.getNotPrimaryAlignmentFlag()){
                result[0]++;

                // check to see if it's aligned:
                if(!read.getReadUnmappedFlag()){
                    result[1]++;
                    if(!read.getDuplicateReadFlag()){
                        result[5]++;
                    }

                    Short alnDistance = read.getShortAttribute("NM");

                    if(alnDistance == null && nmErrorLogCount < 11)
                    {
                        System.out.println("WARNING: no NM attribute found for read " + read);
                        nmErrorLogCount++;
                    }
//
//                  System.out.println("Dist : "+ alnDistance);
                    if (alnDistance == null || alnDistance < DISTANCE_THRESHOLD){
                        result[3]++;
                    }
                }

            }// end of not primary alignment
        }// end of read loop

        referenceSam.close();

        result[2]  = result[1] / result[0];
        result[4]  = result[3] / result[0];

        return result;
    }

    private void pairBWA(String end1SAI, String end2SAI, String end1Fastq, String end2Fastq, String outfile) throws IOException, InterruptedException {
        System.out.println("Running BWA sampe");
        ProcessBuilder pb = new ProcessBuilder(BWA, "sampe",reference, end1SAI, end2SAI, end1Fastq,  end2Fastq);

        System.out.println("Command: " + Tools.collectionToString(pb.command()));
        //pb.redirectErrorStream(true);

        Process pro = pb.start();

        // create thread to read out the standard out
        IOTools.StdOGobbler gobbler = new IOTools.StdOGobbler(pro);
        IOTools.StdErrGobbler errGobbler = new IOTools.StdErrGobbler(pro);

        gobbler.setToFile(outfile);
        gobbler.start();
        errGobbler.start();

        //System.out.println("Waiting on gobbler");
        gobbler.join();
        errGobbler.join();
        //System.out.println("Waiting on process");
        pro.waitFor();
        //System.out.println("Done.");

        System.out.println("Call to BWA complete");
        //System.out.println(IOTools.fileToString(outfile));
    }

    private void singleBWA(String end1SAI, String end1Fastq, String outfile) throws IOException, InterruptedException {
        System.out.println("Running BWA sampe");
        ProcessBuilder pb = new ProcessBuilder(BWA, "samse",reference, end1SAI, end1Fastq);

        System.out.println("Command: " + Tools.collectionToString(pb.command()));
        //pb.redirectErrorStream(true);
        Process pro = pb.start();

        // create thread to read out the standard out
        IOTools.StdOGobbler gobbler = new IOTools.StdOGobbler(pro);
        IOTools.StdErrGobbler errGobbler = new IOTools.StdErrGobbler(pro);
        gobbler.setToFile(outfile);
        gobbler.start();
        errGobbler.start();
        //System.out.println("Waiting on gobbler");
        gobbler.join();
        //System.out.println("Waiting on process");
        errGobbler.join();
        pro.waitFor();
        //System.out.println("Done.");

        System.out.println("Call to BWA complete");
        //System.out.println(IOTools.fileToString(outfile));
    }

    public void getFastQ(String inputBam, String out1, String out2){

        String tmpCommand = "";

        if (new File("/broad/hptmp").exists()){
            System.setProperty("java.io.tmpdir", "/broad/hptmp");

        }

        String[] toFQArgs = null;
        if(out2 == null)
        {
            toFQArgs = new String[]{"I="+inputBam, "F="+out1, "VALIDATION_STRINGENCY=SILENT"};
        }
        else
        {
            toFQArgs = new String[]{"I="+inputBam, "F="+out1, "F2="+out2, "VALIDATION_STRINGENCY=SILENT"};
        }
        SamToFastq stfq = new SamToFastq();

        stfq.VALIDATION_STRINGENCY = SAMFileReader.ValidationStringency.SILENT;
        stfq.instanceMain(toFQArgs);

    }


    private void runBWA(String fastq,  String outfile) throws IOException, InterruptedException {
        System.out.println("Running BWA on " + fastq);
        ProcessBuilder pb = new ProcessBuilder(BWA, "aln",reference, fastq);

        System.out.println("Command: " + Tools.collectionToString(pb.command()));
        //pb.redirectErrorStream(true);
        Process pro = pb.start();

        // create thread to read out the standard out
        IOTools.StdOGobbler gobbler = new IOTools.StdOGobbler(pro);
        IOTools.StdErrGobbler errGobbler = new IOTools.StdErrGobbler(pro);
        gobbler.setToFile(outfile);
        gobbler.start();
        errGobbler.start();

        //System.out.println("Waiting on gobbler");
        gobbler.join();
        errGobbler.join();
        //System.out.println("Waiting on process");
        pro.waitFor();
        //System.out.println("Done.");

        System.out.println("Call to BWA complete");
        //System.out.println(IOTools.fileToString(outfile));
    }

    public int getCount() {
        return (int)result[1];
    }

    public float getRatio(){
        return result[2];
    }

    public int getTotalReads() {
        if (dSampleRate> 0d){
            return (int)(((double)getCount())/ dSampleRate);
        }
        return getCount();
    }

    public int getUniqueReads() {
        if (dSampleRate> 0d){
            return (int)(((double)result[5])/ dSampleRate);
        }

        return (int)result[5];
    }

    /**
     * values of zero or one disable dounsampling (off by default)
     *
     * @param dSampleRate
     */
    public void setdSampleRate(double dSampleRate) {
        this.dSampleRate = dSampleRate;
    }

}
