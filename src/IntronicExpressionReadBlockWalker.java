package org.broadinstitute.cga.rnaseq.gatk;

import net.sf.samtools.SAMRecord;
// import htsjdk.samtools.SAMRecord;
import org.broadinstitute.cga.rnaseq.Transcript;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;

import java.io.*;
import java.util.*;


/**
 * Created by the Broad Institute, Cancer Genome Analysis Group
 * Author: David S. DeLuca
 * User: ddeluca
 * Date: 3/31/11
 * Time: 10:30 AM
 *
 * This adaptation of the CountReadWalker tracks a series of metrics for RNA-seq QC
 *, DataSource.REFERENCE_BASES
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class IntronicExpressionReadBlockWalker extends CountReadBlockMetricsWalker {
    public static final int LOWER_READ_COUNT_CUTOFF =  5;
    
    LinkedList<TranscriptMetrics> transcripts = new LinkedList<TranscriptMetrics>();
    BufferedWriter out;
    BufferedWriter iOut;
    BufferedWriter xOut;
    String lastContig ="";

    String intronReportFile = null;
    String intronOnlyReportFile = null;
    String exonOnlyReportFile = null;
    int totalExonQuantifiedReads = 0;

    @Override
    public void initialize() {
        super.initialize();
        try {  // open outputfile
            this.intronReportFile = super.OUT_FILE+".intronReport.txt";
            out = new BufferedWriter(new FileWriter(this.intronReportFile));
            // Transcript ID | gene name | Number of reads in exons | Total exon length | Exon reads per base | Number of reads in introns |
            // Total intron length | Intron reads per base | Peak intron index | Peak intron read count | Peak intron length | Peak intron reads per base
            // Duplicate exon reads | Split reads
            out.write("Transcript\tGene_Name\tExon_Reads\tExon_Length\tExon_r/b\tIntron_Reads\tIntron_Length\tIntron_r/b\tPeak_Intron_Idx\tPeak_Intron_Reads\tPeak_Intron_Length\tPeak_Intron_r/b\tDuplicate_Exon_Reads\tSplit_Reads\n");

            intronOnlyReportFile = super.OUT_FILE+".intronReport.txt"+"_intronOnly.txt";
            iOut = new BufferedWriter(new FileWriter(intronOnlyReportFile));
            // Intron ID | Transcript ID | Gene name| Number of reads in introns | Total intron length | Intron reads per base |
            // Flanking exons reads | Flanking exons length
            iOut.write("Intron\tTranscript\tGene_Name\tIntron_Reads\tIntron_Length\tIntron_r/b\tFlanking_Exons_Reads\tFlanking_Exons_Length\n");

            exonOnlyReportFile = super.OUT_FILE+".intronReport.txt"+"_exonOnly.txt";
            xOut = new BufferedWriter(new FileWriter(exonOnlyReportFile));
            // Exon ID | Transcript ID | Gene name | Number of reads in exons | Total exon length | Duplicate exon reads | Split reads
            xOut.write("Exon\tTranscript\tGene_Name\tExon_Reads\tExon_Length\tDuplicates\tSplit_Reads\n");
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException("Bad IO for "+ super.OUT_FILE+".intronReport.txt file", e);
        }
    }

    @Override
    protected void makeRefSeqDerviedCounts(SAMRecord read, ArrayList<RefSeqFeature> refSeqs, GenomeLoc readLoc) {
        //check to see if we've switched contigs:
        if (!this.lastContig.equals(read.getReferenceName())) {
            lastContig = read.getReferenceName();
            shedAll();
        }

        boolean wasIntragenicRead = false; // if it's ever in a gene, it's tranagenic
        boolean wasExonicRead = false; // if it's ever an exon, its' an exon
        boolean wasIntronRead = false; // if it's ever an exon, then it's not an intron
        boolean wasJunction = true;   // if it's ever not a junction, then it's not a junction
        boolean transcriptWasPlus = false;
        boolean transcriptWasMinus = false;
        int fragmentLength = 0;
        ArrayList<GenomeLoc> blocks = getBlocksForRead(read);

        // check block for split read:
        boolean splitRead = false;
        int lastBlockEnd = -1;

        for (GenomeLoc blockLoc : blocks) {
            if (lastBlockEnd > 0 && !splitRead) {
                splitRead = ((blockLoc.getStart() - lastBlockEnd) > 100);
            }
            lastBlockEnd = blockLoc.getStop();
        }
        
        for (RefSeqFeature refSeq: refSeqs) {
            int intronNumber = -1;
            boolean thisTranscriptWasExonicRead = false;
            boolean thisTranscriptWasIntronRead = false;

            // sanity check:
            if (!refSeq.overlapsP(readLoc)){
                System.out.println("WHY DOESN'T THIS READ OVERLAP IT'S OWN ROD ENTRY !?!??!!?");
                System.out.println("\tReadLoc:"+readLoc);
                System.out.println("\trefSeqLoc:"+refSeq.getLocation());
                System.out.println("\tRead: " + read);
                System.out.println("\tTransc:" + refSeq.getTranscriptId());
            }

            // calculations for fragmentLength:
            if (fragmentLength == 0 && read.getReadPairedFlag() && read.getProperPairFlag() && thisMismatchCount == 0) {
                fragmentLength = getFragmentLength(refSeq,read,readLoc);
                if (fragmentLength > 0) {
                    super.fragmentLengths.add(fragmentLength);
                }
            }

            HashMap<Integer,Float> splitDosage = null;
            if (splitRead) splitDosage = new HashMap<Integer,Float>();

            for (GenomeLoc blockLoc : blocks) {
        //                        System.out.println("\trefseq:"+refSeq.getGeneName());
        //                        System.out.println("\trefseq:" + refSeq.getTranscriptId());
        //                        System.out.println("\trefseq:"+refSeqLoc);

                //GenomeLoc readLoc = ref.getLocus();
                //intragenic++;// we are in a gene, otherwise, we would be in the ROD
                wasIntragenicRead = true;
                //Performance overlapps = new Performance("Overlap perf: " , Performance.Resolution.milliseconds);
                intronNumber = overlapsIntron(refSeq,blockLoc); // returns intron position if block is on an exon/intron junction

//                System.out.println("\tBLOCK:\t " + blockLoc.toString());
//                System.out.println("\t\tintron:"+ intronNumber);
                if (intronNumber < 0) {
                    thisTranscriptWasExonicRead=true;
                } else {
                    // if we're in a gene but not in an exon ... it must be an intron? it could be UTR :(
                    //intronOrUTR++;
                    thisTranscriptWasIntronRead = true;
        //                            System.out.println("\tNOT Exonic");
                }
                if (splitRead && splitDosage !=null) {
                    if (intronNumber < 0) {
                        Float dosage = splitDosage.get(intronNumber);
                        if(dosage == null) {
                            dosage  = (float)blockLoc.size() / (float)read.getReadLength();
                        } else {
                            dosage  = dosage + (float)blockLoc.size() / (float)read.getReadLength();
                        }
                        splitDosage.put(intronNumber, dosage);
                    } else {
                        splitDosage = null;
                    }
                }
//                if (refSeq.overlapsCodingP(blockLoc)){
//                    //coding++;
//                    wasCodingRead = true;
//                            System.out.println("\tIs Coding");
//                }
                //overlapps.outputIfTookTime();
            }// end of blocks

            // to do with this RefSeq entry:
            try {
                // if a readblock had overlapped an exon/intron junction, it was assigned to the intron
                // (both as intron number as well as transcriptWasExonicRead)
                // if the read consists of multiple blocks, then the exon/intron chosen is done so arbitrarly: the latter exon will be chosen, because the loop above just erases the previous exon number
                // this means that if an exon is not included in this isoform, and the read spans those exons
                // it will appear here as an intronic read and be counted as an intron for the individual
                // intron metrics. Futhermore, if a read lands completely in the exon from a different isoform
                // the same issue will happen - it will be counted as an intronic read.
                // HOWEVER, NOTE: that this issue will affect the transcript-specific counts that are saved in the
                // exons file, but not the global counts reported in the final metrics.
                //todo probably, it's best just to run this on a gene model annotation and not a transcript model annotation

                boolean keep = true;
                if (super.strictExonReads) {
                    if(read.getReadPairedFlag()) {
                        Short alnDistance = read.getShortAttribute("NM");
                        keep = read.getProperPairFlag() && alnDistance!=null && alnDistance <=6 && read.getMappingQuality() >3;
                    }
                }
                if (keep) {
                    addRefSeqToMetrics(refSeq, intronNumber, read.getDuplicateReadFlag(),splitRead, splitDosage);
                }
            } catch (RuntimeException e) {
                System.out.println("Exception adding metrics to linked list");
                System.out.println("RefSeq: " + refSeq.toString());
                System.out.println("Read: " + read.toString());
                System.out.println("Intron number: " + intronNumber);
                throw e;
            }

            //JUNCTION:  if it is ever purely exonic, then it was not a junction
            if (thisTranscriptWasExonicRead && !thisTranscriptWasIntronRead){
//                System.out.println("\t\tBlock was junction");
                wasJunction = false;
            } else {
//                System.out.println("\t\tBlock was NOT a junction");
            }
            if (thisTranscriptWasExonicRead) wasExonicRead = true;
            if (thisTranscriptWasIntronRead) wasIntronRead = true;

            shedMetrics(read.getAlignmentStart());

            //check for strand specificity:
            if (refSeq.getStrand() == 1) {  // end1 reads is sense strand
                transcriptWasPlus = true;
            } else {
                transcriptWasMinus = true;
            }
        }// end of RefSeq entry iteration
        if (wasIntragenicRead) {
            intragenic++;
            if (wasExonicRead && !wasJunction){
                exonic++; // note in this case, we're defining junction reads as intronic
                if (splitRead) {
                    this.splitReads++;
                }
            } else {
                intronOrUTR++; // was an intron or junction
            }
        }

        // count strand specificity metrics:
        super.countStrandSpecificMetrics(read, transcriptWasMinus,transcriptWasPlus);
        //if (wasCodingRead) coding++;
        // junction reads? in this case we're assuming  that they're intronic
    }

    private int getFragmentLength(RefSeqFeature refSeq, SAMRecord read, GenomeLoc readLoc) {
        int exon = inExon(refSeq,readLoc);
        if (exon > 0) {
            if (!read.getMateUnmappedFlag()) {
                // if the mate is mapped:
                if (read.getReferenceIndex().equals(read.getMateReferenceIndex())) {
                    // we're on the same chromosome
                    int mateExon = inExon(refSeq,read.getMateAlignmentStart());
                    if (mateExon == exon) {
                        // we can use these reads to figure out the fragment size
                        if (read.getFirstOfPairFlag()) {
                            return (read.getMateAlignmentStart() + read.getReadLength()) - readLoc.getStart();
                        } else {
                            return readLoc.getStop() - read.getAlignmentStart();
                        }
                    }
                }
            }
        }
        return 0;
    }

    private void shedAll() {
        for (TranscriptMetrics met:  transcripts){
            outputMetric(met);
            if (met.totalExonReadCount >= LOWER_READ_COUNT_CUTOFF) {
                geneLog.add(met.name);
                transcriptLog.add(met.transcriptId);
            }
        }
        transcripts.clear();
    }

    @Override
    public void onTraversalDone(Integer result) {
        super.onTraversalDone(result);

        try {
            shedAll();// output any remaining transcripts
            out.close();
            System.out.println("Finished writing " + this.intronReportFile);

            //make a new gct file for RPKM: fragments per kilobase of exon per million fragments mapped
            // first, we have to load all the data from the intronreport file and calc RPKMs
            BufferedReader in =  new BufferedReader (new FileReader(this.intronReportFile));

            String line = in.readLine(); // skip header
            line = in.readLine();

            HashMap<String,String> transcriptRPKMMap = new HashMap<String,String>();
            while (line != null) {
                String [] split = line.split("\\t");
                String id = split[0];
                String name = split[1];
                int exonReads= Integer.parseInt(split[2]);
                int exonLength = Integer.parseInt(split[3]);

                float rpkm = Transcript.getRPKM(exonReads, exonLength, this.totalExonQuantifiedReads); // strict mode support

                transcriptRPKMMap.put(id,id+'\t'+name+'\t'+rpkm);
                line = in.readLine();
            }
            in.close();

            // read in the orginal refseq file so that we can output all the values in this original order
            // making sure that we also output zeros if the transcript was never seen by the walker
            BufferedReader refSeqIn = new BufferedReader (new FileReader(super.RefseqFileName));
            line = refSeqIn.readLine();
            ArrayList<String> dataLines = new ArrayList<String>();
            while (line != null){
                String [] split = line.split("\\t");
                String id = split[1];
                String name = split[12];
                String dataLine = transcriptRPKMMap.get(id);
                if (dataLine == null){
                    dataLine = id+'\t'+name+"\t0"; // value is zero
                }
                dataLines.add(dataLine);
                line = refSeqIn.readLine();
            }
            refSeqIn.close();

            // now we can output the RPKM dat in GCT format:
            out = new BufferedWriter (new FileWriter(super.OUT_FILE + ".rpkm.gct"));
            out.write("#1.2\n");
            out.write("" + dataLines.size() + "\t1\n"); // just one sample
            out.write("NAME\tDescription\tSAMPLE\n");
            // so that no transcripts get lost
            for (String data: dataLines){
                out.write(data);
                out.write('\n');
            }
            out.close();

            // introns.  first we close the file that contains the intron data
            iOut.close();

            System.out.println("Finished writing " + this.intronOnlyReportFile+", now creating RPKM values for introns ..");

            // INTRON ONLT RPKMS
            in = new BufferedReader (new FileReader(this.intronOnlyReportFile));

            line = in.readLine(); // skip header
            line = in.readLine();

            StringBuilder dataSection = new StringBuilder();

            //todo make this file load the RefseqFileName and conform to its order
            int dataCount = 0;
            while (line != null){ // this loop reads in all the intron-based info from the intronOnly file and then creates a large data string
                String [] split = line.split("\\t"); // has form: Intron	Transcript	Gene_Name	Intron_Reads	Intron_Length	Intron_r/b
                String id = split[0];
                String name = split[1] +";" + split[2]; // in GCT file this is actually the description field
                int intronReads= Integer.parseInt(split[3]);
                int intronLength = Integer.parseInt(split[4]);
                if (intronLength>0){
                    float rpkm = Transcript.getRPKM(intronReads, intronLength, this.totalExonQuantifiedReads); //strictmode support
                    dataSection.append(id).append('\t').append(name).append('\t');
                    dataSection.append(String.format("%.3f", rpkm)).append('\n');
                    dataCount++;
                }
                line = in.readLine();
            }
            in.close();

            // now we can output the RPKM dat in GCT format:
            out = new BufferedWriter (new FileWriter(super.OUT_FILE + ".introns.rpkm.gct"));
            out.write("#1.2\n");
            out.write("" + dataCount + "\t1\n"); // just one sample
            out.write("NAME\tDescription\tSAMPLE\n");
            out.write(dataSection.toString());
            out.close();

            xOut.close();

            FileWriter fout = new FileWriter(super.OUT_FILE + ".totalExonQuantifiedReads");
            fout.write(""+totalExonQuantifiedReads);
            fout.close();
            //todo reformat exons into GCT
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException("Could not close file "+ super.OUT_FILE+".intronReport.txt ", e);
        }
    }


    private void shedMetrics(int readStart) {
        for (Iterator<TranscriptMetrics> it = transcripts.iterator(); it.hasNext();){
            TranscriptMetrics met = it.next();
            if (met.getEnd() < readStart){
                outputMetric(met);
                if (met.totalExonReadCount >= LOWER_READ_COUNT_CUTOFF){
                    geneLog.add(met.name);
                    transcriptLog.add(met.transcriptId);
                }
                it.remove();
            }
        }
    }

    /**
     * The output format is tab delimited with fields:
     *    Transcript ID
     *    Number of reads in exons
     *    Total exon length
     *    Exon Reads per base
     *    Number of reads in introns
     *    Total intron length
     *    Intron Reads per Base
     *    Peak intron index
     *    Peak intron read count
     *    Peak intron length
     *    Peak intron reads per base
     *
     *    Total Duplicate Exon Reads
     *
     *
     * @param met
     */
    private void outputMetric(TranscriptMetrics met)  {
    //"Transcript\tGene_Name\tExon_Reads\tExon_Length\tExon_r/b\tIntron_Reads\tIntron_Length\tIntron_r/b\tPeak_Intron_Idx\tPeak_Intron_Reads\tPeak_Intron_Length\tPeak_Intron_r/b\n")
        try{
            out.write(met.transcriptId); out.write('\t');
            out.write(met.name); out.write('\t');
            out.write(""+met.totalExonReadCount); out.write('\t');
            out.write(""+met.totalExonLength); out.write('\t');
            out.write(String.format("%.3f",met.getExonReadsPerBase())); out.write('\t');
            out.write(""+met.totalIntronReadCount); out.write('\t');
            out.write(""+met.totalIntronLength); out.write('\t');
            out.write(String.format("%.3f",met.getIntronReadsPerBase())); out.write('\t');
            out.write(""+met.getPeakIndex()); out.write('\t');
            out.write(""+met.getPeakReadCount()); out.write('\t');
            out.write(""+met.getPeakLength()); out.write('\t');
            out.write(String.format("%.3f",met.getPeakReadsPerBase())); out.write('\t');
            out.write(""+met.totalDuplicateExonReadCount); out.write('\t');
            out.write(""+met.totalSplitReadCount); out.write('\n');
            
            // intron only: Intron	Transcript	Gene_Name	Intron_Reads	Intron_Length	Intron_r/b
            for (int i = 0 ; i < met.intronReadCount.length; i++){
                iOut.write(met.transcriptId+"_"+i); iOut.write('\t');
                iOut.write(met.transcriptId); iOut.write('\t');
                iOut.write(met.name); iOut.write('\t');
                iOut.write(""+met.intronReadCount[i]); iOut.write('\t');
                iOut.write(""+met.intronLengths[i]); iOut.write('\t');
                iOut.write(String.format("%.3f", ((float)met.intronReadCount[i] /(float)met.intronLengths[i]))); iOut.write('\t');
                // flanking exon reads + flanking exon lengths

                float exonReads = 0; int exonLength = 0;
                if (i >0){
                    // left exon:
                    exonReads= met.exonReadCount[i-1];
                    exonLength = met.exonLengths[i-1];
                }
                if (i < met.exonReadCount.length){
                    // right exon
                    exonReads += met.exonReadCount[i];
                    exonLength += met.exonLengths[i];
                }
                iOut.write(""+exonReads+"\t"+exonLength);
                iOut.write('\n');
            }

            // exon only: "Exon\tTranscript\tGene_Name\tExon_Reads\tExon_Length\tDuplicates\n"
            for (int i = 0 ; i < met.exonReadCount.length; i++){
                xOut.write(met.transcriptId+"_"+i); xOut.write('\t');
                xOut.write(met.transcriptId); xOut.write('\t');
                xOut.write(met.name); xOut.write('\t');
                xOut.write(String.format("%.3f",met.exonReadCount[i])); xOut.write('\t');
                xOut.write(""+met.exonLengths[i]); xOut.write('\t');
                xOut.write(""+met.exonDuplicateReadCount[i]); xOut.write('\t');
                xOut.write(""+met.exonSplitReadCount[i]);
                xOut.write('\n');
            }
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException("Could not close file "+ super.OUT_FILE+".intronReport.txt ", e);
        }
    }

    /**
     * This method determines which exon or intron the block (portion of a read) is in.
     * If it crosses an exon/intron junction, then it is called as an intronic read
     *
     *
     * @param refSeq
     * @param blockLoc
     * @return -1 if it's in an exon 1, -2 for exon 2, etc
     *           otherwise it returns the number of the intron that this read overlaps with
     *           Note: a value of zero refers to the space before the first exon
     */
    private int overlapsIntron(RefSeqFeature refSeq, GenomeLoc blockLoc) {
        int i = 0;
        for (GenomeLoc exon: refSeq.getExons()){
            i++;
            if (blockLoc.getStart() < exon.getStart()){
                // our read starts before this exon, so return (could be a junction position)
                return i-1;
            } else if (blockLoc.getStop() <= exon.getStop()) {  // read is in the exon
                return 0-i; // actual exon number, but as a negative
            } else if (blockLoc.getStart() <= exon.getStop()) {
                // this read is at the junction at the end of this exon
                return i;
            }
            //at this point this read has to either be in the next intron or down past the next exon
            // so we'll pick  it up in a subsequent iteration of this loop
        }
        return i; // this is beyond the last exons;
    }


    private void addRefSeqToMetrics(RefSeqFeature refSeq, int intronNumber, boolean wasDuplicate, boolean wasSplit, HashMap<Integer,Float> dosage) {
        String id = refSeq.getTranscriptId();

        for (TranscriptMetrics tm: this.transcripts){
            if (tm.transcriptId.equals(id)){
                // add to this intron.
                tm.countIntronExon(intronNumber, wasDuplicate, wasSplit, dosage);
                return;
            }
        }
        // no transcript found, add a new one:
        TranscriptMetrics newMetrics = new TranscriptMetrics(refSeq);
        newMetrics.countIntronExon(intronNumber,wasDuplicate, wasSplit, dosage);
        this.transcripts.add(newMetrics);
    }


    private int inExon(RefSeqFeature refSeq, GenomeLoc readLoc) {
        int i = 0;
        for (GenomeLoc exon: refSeq.getExons()) {
            i++;
            if (readLoc.getStart() >= exon.getStart() && readLoc.getStop() <= exon.getStop()) {
                return i;
            }
        }
        return 0; // this is beyond the last exons;
    }


    private int inExon(RefSeqFeature refSeq, int pos) {
        int i = 0;
        for (GenomeLoc exon: refSeq.getExons()){
            i++;
            if (pos >= exon.getStart() && pos <= exon.getStop()) {
                return i;
            }
        }
        return 0; // this is beyond the last exons;
    }


    class TranscriptMetrics {
        int totalExonLength;   // constructor
        int totalIntronLength; // constructor
        //int totalLength;

        int totalExonReadCount =0; //incremented
        int totalIntronReadCount=0;//incremented

        int [] intronReadCount; // filled with zeros in construtor
        float [] exonReadCount; // filled with zeros in construtor

        int totalDuplicateExonReadCount =0; //incremented
        int totalSplitReadCount=0;
        int [] exonDuplicateReadCount; // filled with zeros in construtor
        int [] exonSplitReadCount; // filled with zeros in construtor
        //int [] exonReadCount;

        String transcriptId, name; // constructor
        private int end;     // constructor
        private int peakIndex = 0; // incremented
        private int[] intronLengths; // constructor
        private int[] exonLengths; // constructor

        TranscriptMetrics (RefSeqFeature refSeq){
            // TRANSCRIPT ID
            this.transcriptId = refSeq.getTranscriptId();
            this.name = refSeq.getGeneName();

            // intron read counts initialzied to zero
            this.intronReadCount = new int[refSeq.getExons().size()+1];
            Arrays.fill(this.intronReadCount,0);
            this.exonReadCount = new float[refSeq.getExons().size()];
            Arrays.fill(this.exonReadCount,0);

            this.exonDuplicateReadCount = new int[refSeq.getExons().size()];
            Arrays.fill(this.exonDuplicateReadCount,0);

            this.exonSplitReadCount = new int[refSeq.getExons().size()];
            Arrays.fill(this.exonSplitReadCount,0);

            // intron lengths calculated from exon lengths:
            this.intronLengths = new int[refSeq.getExons().size()+1];
            this.exonLengths = new int[refSeq.getExons().size()];
            
            int i = 0;
            int lastIntronStart = refSeq.getStart();
            totalExonLength = 0;
            totalIntronLength =0;
            for (GenomeLoc exon: refSeq.getExons()){
                int intronLength = exon.getStart() - lastIntronStart;
                if (intronLength < 0){
                    System.out.println("Warning, bad boundary: " + transcriptId);
                    intronLength = 0;
                }
                this.intronLengths[i] = intronLength;
                totalIntronLength+=this.intronLengths[i];

                this.exonLengths[i] = exon.getStop() - exon.getStart() + 1;
                totalExonLength+= exonLengths[i];
                lastIntronStart = exon.getStop()+1;
                i++;
            }
            int lastLength = refSeq.getEnd() - lastIntronStart +1;
            if (lastLength < 0 ){
                System.out.println("Warning bad boundary (at end) " +transcriptId);
                lastLength =0;
            }
            this.intronLengths[i] = lastLength;
            totalIntronLength+=intronLengths[i];

            // refseq end:
            this.end = refSeq.getEnd();
        }


        public void countIntronExon(int intronNumber, boolean wasDuplicate, boolean wasSplit, HashMap<Integer,Float> dosage) {
            if (intronNumber < 0) {  // this is an exon
                totalExonQuantifiedReads++;
                totalExonReadCount++;
                if (wasSplit && dosage!=null){  // using dosage
                    for (Integer dosageIntrNum: dosage.keySet()) {
                        exonReadCount[-1*dosageIntrNum-1]+= dosage.get(dosageIntrNum);
                    }
                } else {
                    exonReadCount[-1*intronNumber-1]++;
                }
                if (wasDuplicate) {
                    totalDuplicateExonReadCount++;
                    exonDuplicateReadCount[-1*intronNumber-1]++;
                }
                if (wasSplit) {
                    totalSplitReadCount++;
                    exonSplitReadCount[-1*intronNumber-1]++;
                }
            } else {
                intronReadCount[intronNumber]++;
                totalIntronReadCount++;
                if (intronReadCount[intronNumber] > intronReadCount[peakIndex]){
                    peakIndex = intronNumber;
                }
            }
        }

        public int getEnd() {
            return end;
        }

        public float getExonReadsPerBase() {
            return (float)this.totalExonReadCount / (float) this.totalExonLength;
        }

        public float getIntronReadsPerBase() {
            if (this.totalIntronLength > 0)
                return (float)this.totalIntronReadCount / (float) this.totalIntronLength;
            else
                return -1f;
        }

        public float getPeakReadsPerBase() {
            if (intronLengths[peakIndex] > 0)
                return (float) intronReadCount[peakIndex] / (float) this.intronLengths[peakIndex];
            else
                return -1f;
        }

        public int getPeakIndex() {
            return peakIndex;
        }

        public int getPeakReadCount() {
            return intronReadCount[peakIndex];
        }

        public int getPeakLength() {
            return this.intronLengths[peakIndex];
        }
    }
}
