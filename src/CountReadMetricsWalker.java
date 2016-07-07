package org.broadinstitute.cga.rnaseq.gatk;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.cga.tools.ObjectCounter;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.SeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqCodec;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

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
 *
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CountReadMetricsWalker extends ReadWalker<Integer,Integer> {
    @Argument(fullName = "outfile_metrics", shortName = "OM", doc="The destination file for the metrics", required=true)
    protected String OUT_FILE = null;

    @Argument(fullName="refseq", shortName="refseq", doc="Name of RefSeq transcript annotation file.", required=false)
    protected String RefseqFileName = null;

    @Argument(fullName = "strict_exon_reads", shortName = "strict", doc="Apply filters to reads before quantifying", required=true)
    protected boolean strictExonReads = false;

    private LocationAwareSeekableRODIterator refseqIterator=null;

    public static final int CHIMERA_DISTANCE = 2000000; // 2 mb

    GenomeLoc lastInterval = null;
    GenomeLoc longestInterval = null;
    ArrayList<RefSeqFeature> lastAnnotationList = null;
    ArrayList<RefSeqFeature> longestRefSeqList = null;

    int totalReads = 0;
    int mappedReads = 0;
    int altAlignments = 0;
    int failed = 0;
    int readLength = 0;

    protected int mappedUniqueReads = 0;
    int duplicates = 0;
    //int unique = 0;
    int mappedDuplicates;
    int overlaps = 0;
    //int coding = 0;
    int intragenic = 0;
    int exonic = 0;
    int splitReads= 0;
    int intergenic = 0;
    int intronOrUTR = 0;
    int droppedByLength = 0;
    int end1Sense = 0;
    int end1Antisense = 0;
    int end2Sense = 0;
    int end2Antisense = 0;

    int end1Mapped = 0;
    int end2Mapped = 0;

    int thisMismatchCount = 0;
    ObjectCounter<Integer> fragmentLengths = new ObjectCounter<Integer>();
    HashSet<String> transcriptLog= new HashSet<String>();
    HashSet<String> geneLog = new HashSet<String>();

    final int MAX_READ_LEGNTH = 100000;

    int longestReadLength = 0;
    int maxFeatureSize = 0;

    long totalBases = 0;
    long mismatchBases = 0;

    long end1MismatchBases = 0;
    long end2MismatchBases = 0;

    long end1TotalBases = 0;
    long end2TotalBases = 0;

    int totalMappedPairs = 0;
    int chimericPairs = 0;
    int unpairedReads = 0;

    BufferedWriter chOut = null;
/*
public List<Object> getReferenceMetaData(final String name) {
       RODRecordList list = getTrackDataByName(name, true);
       List<Object> objects = new ArrayList<Object>();
       if (list == null) return objects;
       for (GATKFeature feature : list)
           objects.add(feature.getUnderlyingObject());
       return objects;
   }
*/

    @Override
    public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
        if (read.getNotPrimaryAlignmentFlag()) {
            altAlignments++;
        } else {
            if (read.getReadFailsVendorQualityCheckFlag()) {
                failed ++;
            }
        }

        //Short alnDistance = read.getShortAttribute("NM");
        //todo: if (!read.getNotPrimaryAlignmentFlag() && read.getProperPairFlag() && alnDistance!=null && alnDistance <=6 && read.getMappingQuality() >0){

        if (!read.getNotPrimaryAlignmentFlag() && !read.getReadFailsVendorQualityCheckFlag()){
            totalReads++;
            if(!read.getReadPairedFlag()){
                unpairedReads++;
            }
            if (read.getReadLength() > readLength){
                readLength = read.getReadLength();
            }
            if (read.getDuplicateReadFlag()){
                duplicates++;
            }else{
                //unique++; // removed because it is not true that it is unique, since "not duplicate" set also contains unaligned reads
            }
            if (!read.getReadUnmappedFlag()){
                mappedReads++;
                if (read.getReadPairedFlag()){
                    if (read.getFirstOfPairFlag()){
                        end1Mapped++;
                    }else{
                        end2Mapped++;
                    }
                }
                // count number of mismatched bases:
                countMismatchedBases(ref, read);
                if (read.getDuplicateReadFlag()){
                    mappedDuplicates++;
                }else{
                    mappedUniqueReads++;
                }
                //Performance perf = new Performance("Total Map Call: ", Performance.Resolution.milliseconds);
                GenomeLoc readLoc =  getToolkit().getGenomeLocParser().createGenomeLoc(read); // ref.getLocus();
                int length = readLoc.getStop() - readLoc.getStart() +1;

                if (length > MAX_READ_LEGNTH){
                    droppedByLength++;
                    return 1;
                }
                if (length > longestReadLength){
                    longestReadLength = length;
                }
//                System.out.println("Read:\t"+read.toString());
//                System.out.println(" at location:\t"+location);
//                System.out.println(" is duplicate:\t" + read.getDuplicateReadFlag());
                ArrayList<RefSeqFeature> refSeqs =  getRodList(readLoc);

                //check chimeric pairs
                if (read.getReadPairedFlag() && !read.getMateUnmappedFlag()){  // if mate is also mapped
                    if (read.getFirstOfPairFlag()) {
                        totalMappedPairs++;
                    }
                    if (read.getReferenceIndex() != read.getMateReferenceIndex() ||
                            (Math.abs(read.getAlignmentStart()-read.getMateAlignmentStart()) > CHIMERA_DISTANCE)) {
                        chimericPairs++;
                        String transcript = "";
                        String gene = "";
                        if (refSeqs != null && refSeqs.size() > 0){
                            RefSeqFeature r = refSeqs.get(0);
                            transcript = r.getTranscriptId();
                            gene = r.getGeneName();
                        }
                        try {  // write chimeric pairs
                            chOut.write(transcript+"\t"+gene+"\t"+readLoc.toString()+"\t"+read.getMateReferenceName()+":"+read.getMateAlignmentStart()
                                +"\t"+(read.getReferenceIndex() == read.getMateReferenceIndex())+"\n");
                        } catch(IOException e) {
                            throw new RuntimeException(e.getMessage());
                        }
                    }
                }
                
                // sanity check:
                if (read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START){
//                    System.out.println("HEY, READ IS MAPPED BUT HAS NO ALIGNMENT START!!!!!!! " +read.toString());
                }

                // read is mapped, we can now have our way with it
                //List<GATKFeature> features = tracker.getAllCoveringRods();

                if (refSeqs ==null){
                    intergenic++; // no annotation, so we're inbetween genes
//                    System.out.println("\tIntergenic");
                }else{
                    int numOfRefSeqs = refSeqs.size();
                    if (numOfRefSeqs > this.maxFeatureSize){
                        //System.out.println("Largest Num of features:\t"+numOfRefSeqs);
                        this.maxFeatureSize = numOfRefSeqs;
                    }
                    //Collection<RefSeqFeature> refSeqs = this.getRefSeqs(annotationList); // this.condenseRefSeqFeatures(tracker.getAllCoveringRods());
                    //Collection<RefSeqFeature> refSeqs = annotationList;
                    if (numOfRefSeqs == 0){
                        // no refseq track means we're in an intragenic region
                        // could this really happeN?
                        //System.out.println("Could this happen? refSeqs collection with size zero?");
//                        System.out.println("\tItergenic, but within span of previous read");
                        intergenic++;
                    }

                    if (numOfRefSeqs> 1){
//                        System.out.println("\toverlapping transcripts!");
                        overlaps++;
                    }

                    makeRefSeqDerviedCounts(read, refSeqs, readLoc);
                    //forRefs.outputIfTookTime();
                }
                //perf.outputIfTookTime();
            } // end of if alligned
        }// end of not primary alignment check
        return 1;
    }


    private void countMismatchedBases(ReferenceContext ref, SAMRecord read) {
        if (ref == null) {
            System.out.println("Why is my ReferenceContext null?");
            System.out.println("Read length = " + read.getReadLength());
            System.out.println(read.getReadString());
            System.out.println(read.toString());
            System.out.println(read.getCigarString());
            thisMismatchCount = -1;
        }else{
            thisMismatchCount = AlignmentUtils.getMismatchCount(read,ref.getBases(),0).numMismatches;
            mismatchBases+= thisMismatchCount;
            totalBases+=read.getReadLength();
            if (read.getReadPairedFlag()){
                if(read.getFirstOfPairFlag()){
                    end1MismatchBases +=thisMismatchCount;
                    end1TotalBases += read.getReadLength();
                }else{
                    end2MismatchBases +=thisMismatchCount;
                    end2TotalBases += read.getReadLength();
                }
            }
        }
    }


    /**
     * Only called if refSeqs has at least one elelement
     *
     * @param read
     * @param refSeqs
     * @param readLoc
     */
    protected void makeRefSeqDerviedCounts(SAMRecord read, ArrayList<RefSeqFeature> refSeqs, GenomeLoc readLoc) {
        //Performance forRefs = new Performance("RefSeq get location loop",Performance.Resolution.milliseconds);
        boolean wasIntragenic = false;
        boolean wasExonic = false;
        boolean wasIntron = false;
        //boolean wasCoding = false;
        boolean transcriptWasPlus = false;
        boolean transcriptWasMinus = false;

        for (RefSeqFeature refSeq: refSeqs){
            GenomeLoc refSeqLoc = refSeq.getLocation();
    //                        System.out.println("\trefseq:"+refSeq.getGeneName());
    //                        System.out.println("\trefseq:" + refSeq.getTranscriptId());
    //                        System.out.println("\trefseq:"+refSeqLoc);

            //GenomeLoc readLoc = ref.getLocus();
            // sanity check:
            if (!refSeq.overlapsP(readLoc)){
                System.out.println("WHY DOESN'T THIS READ OVERLAP IT'S OWN ROD ENTRY !?!??!!?");
                System.out.println("\tReadLoc:"+readLoc);
                System.out.println("\trefSeqLoc:"+refSeqLoc);
                System.out.println("\tRead: " + read);
                System.out.println("\tTransc:" + refSeq.getTranscriptId());
            }
            //intragenic++;// we are in a gene, otherwise, we would be in the ROD
            wasIntragenic = true;
    //                        System.out.println("\tIs Intragenic ...");
            //Performance overlapps = new Performance("Overlap perf: " , Performance.Resolution.milliseconds);

            if (refSeq.overlapsExonP(readLoc)){
    //                            System.out.println("\tIs Exonic");
                //exonic++;
                wasExonic=true;
            }else{
                // if we're in a gene but not in an exon ... it must be an intron? it could be UTR :(
                //intronOrUTR++;
                wasIntron = true;
    //                            System.out.println("\tNOT Exonic");
            }

//            if (refSeq.overlapsCodingP(readLoc)){
//                //coding++;
//                wasCoding = true;
//    //                            System.out.println("\tIs Coding");
//            }
            //overlapps.outputIfTookTime();
            //sheck for strand specificity:
            if (refSeq.getStrand() == 1){
                // end1 reads is sense strand
                transcriptWasPlus = true;
            }else{
                transcriptWasMinus =true;
            }
        } // end of iterating overa all transcript models for this read

        if (wasIntragenic) intragenic++;
        if (wasExonic) exonic++;
        if (wasIntron) intronOrUTR++;
//        if (wasCoding) coding++;
        countStrandSpecificMetrics(read, transcriptWasMinus,transcriptWasPlus);
    }

    protected void countStrandSpecificMetrics(SAMRecord read, boolean transcriptWasMinus, boolean transcriptWasPlus) {
        if (transcriptWasPlus && transcriptWasMinus){
            // hmm. ... what can we do if there were both situations? nothing?
        }else if (!transcriptWasPlus && !transcriptWasMinus){
            // transcript is neither sense nor antisense ....
        }else{
            if (read.getReadPairedFlag()){
                if (read.getFirstOfPairFlag()){ // read is END1
                    if (read.getReadNegativeStrandFlag()){
                        // if the read is on the negative strand, then it is sense if the transcript is annotated to negative strand
                        if (transcriptWasMinus){
                            this.end1Sense++;
                        }else{ // transcript was plus
                            // the read is on the negative strand and transcript is on plus strand then this was an antisense read
                            this.end1Antisense++;
                        }
                    }else{ // read is end1 and on the plus strand
                        if (transcriptWasPlus){
                            this.end1Sense++;
                        }else{
                            this.end1Antisense++;
                        }

                    }
                }
                if(read.getSecondOfPairFlag()){ // this read is End2
                    if (read.getReadNegativeStrandFlag()){
                        // if the read is on the negative strand, then it is sense if the transcript is annotated to negative strand
                        if (transcriptWasMinus){
                            this.end2Sense++;
                        }else{ // transcript was plus
                            // the read is on the negative strand and transcript is on plus strand then this was an antisense read
                            this.end2Antisense++;
                        }
                    }else{ // read is end1 and on the plus strand
                        if (transcriptWasPlus){
                            this.end2Sense++;
                        }else{
                            this.end2Antisense++;
                        }
                    }
                }
            }else{
                //what do we do if these are not paired end reads ...
            }
        }
    }


    /**
     * This method is necessary to overcome a difficulty with the refseqIterator.
     * The problem is, that the interator never wants to encounter a read with an earlier positoin than the last read.
     * This can happen, however, because BAMs are sorted by start position, with no regard to the stop position.
     *
     * @param loc
     * @return
     */
    protected ArrayList<RefSeqFeature> getRodList(GenomeLoc loc) {
        if (lastInterval == null || !lastInterval.getContig().equals(loc.getContig())){
            // first use of interval, so set it and return:
            lastInterval = loc;
            lastAnnotationList = this.getUnderlyingRefSeqs(refseqIterator.seekForward(loc));
            longestInterval = loc;
            this.longestRefSeqList = lastAnnotationList;
            return lastAnnotationList;
        }

        if (loc.equals(lastInterval)){
            return lastAnnotationList;
        }


//        if (loc.getStart() > longestInterval.getStop()){
//            // we're someplace new, so we can rely on the seakForward
//            lastInterval = loc;
//            lastAnnotationList = this.getUnderlyingRefSeqs(refseqIterator.seekForward(loc));
//            longestInterval = loc;
//            longestRefSeqList = lastAnnotationList;
//            return lastAnnotationList;
//        }

        if (loc.getStop() >= longestInterval.getStop()){
            //we've pushed forward, or possibly the start has moved forward, but the end is still the same
            lastInterval = loc;
            lastAnnotationList = this.getUnderlyingRefSeqs(refseqIterator.seekForward(loc));
            this.longestRefSeqList = lastAnnotationList;
            longestInterval = loc;
            return lastAnnotationList;
        }

        // if we make it this far, then the stop position is less than before
        lastInterval = loc; // we have to reclaim the longest interval RODs and then reduce them if some don't overlap
        //lastAnnotationList = refseqIterator.seekForward(longestInterval);
        if (lastAnnotationList == null) return null;
        int wasCount = 0;

        //RefSeqFeature removed = null; // this obj is just for debugging
        this.lastAnnotationList = new ArrayList<RefSeqFeature>();
        for (RefSeqFeature f: this.longestRefSeqList){
            wasCount++;
            if (f.getLocation().overlapsP(loc)){
                this.lastAnnotationList.add(f);
            }else{
                //removed = f;
            }
        }
        /*
        //todo check whether this is correct
        if (wasCount > 0 && lastAnnotationList.size() == 0){
            System.out.println("Reduced annoation list from " + wasCount + " to zero");
            System.out.println("\tThis Loc:\t"+loc);
            System.out.println("\tLongest Loc\t"+longestInterval);
            System.out.println("\tRemoved Trans\t" + removed.getLocation());
        }*/
        return lastAnnotationList;
    }

    private ArrayList<RefSeqFeature> getUnderlyingRefSeqs(RODRecordList annotationList) {
        if (annotationList == null) return null;
        ArrayList<RefSeqFeature> list = new ArrayList<RefSeqFeature>(annotationList.size());
        for (GATKFeature f: annotationList){
            RefSeqFeature rs = (RefSeqFeature)f.getUnderlyingObject();
            list.add(rs);
        }
        return list;
    }


    @Override
    public void initialize() {
        super.initialize();
        MathUtils.ratio(3, 4);

        if ( RefseqFileName != null ) {
            logger.info("Using RefSeq annotations from "+RefseqFileName);

            RMDTrackBuilder builder = new RMDTrackBuilder(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),
                    getToolkit().getGenomeLocParser(),
                    getToolkit().getArguments().unsafe);
            RMDTrack refseq = builder.createInstanceOfTrack(RefSeqCodec.class,new File(RefseqFileName));

            refseqIterator = new SeekableRODIterator(refseq.getHeader(),
                    refseq.getSequenceDictionary(),
                    getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),
                    getToolkit().getGenomeLocParser(),
                    refseq.getIterator());
        }

        if (refseqIterator == null) {
            logger.info("No gene annotations available");
        }
        try {
            chOut = new BufferedWriter(new FileWriter(OUT_FILE+ ".chimericPairs.txt"));
            chOut.write("Transcript Id\tGene Name\tRead Location\tMate Location\tSame Chromosome\n");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public Integer reduceInit() {
        return 0;
    }


    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    @Override
    public void onTraversalDone(Integer result) {
        super.onTraversalDone(result);

        PrintWriter out;
        boolean fileOk = true; // we're going to do some weird handling to output to file if it's there. if
        // no file is given, we'll go to stdout or if the file fails to open
        if (OUT_FILE == null){
            out  = new PrintWriter(System.out);
        } else {
            try {
                out = new PrintWriter (new FileWriter(OUT_FILE+".metrics.tmp.txt"));  // assume that OUT_FILE is <out_dir>/<sample_id>/<sample_id>
            } catch (IOException e){
                e.printStackTrace();
                System.out.println("IOException for file: " + OUT_FILE);
                fileOk=false;
                out = new PrintWriter(System.out);
            }
        }

        /*
        out.println("Total\tMapped\tDuplicates\tOverlaps\tCoding\tIntragenic\tExonic\tIntergenic\tIntronOrUTR\tDropped");

        out.print(totalReads);out.print('\t');
        out.print(mappedReads);out.print('\t');
        out.print(duplicates);out.print('\t');
        out.print(overlaps);out.print('\t');
        out.print(coding);out.print('\t');
        out.print(intragenic);out.print('\t');
        out.print(exonic);out.print('\t');
        out.print(intergenic);out.print('\t');
        out.print(intronOrUTR);out.print('\t');
        out.print(droppedByLength);
*/
        out.println("Total Purity Filtered Reads Sequenced\tAlternative Aligments\tFailed Vendor QC Check\tRead Length"); // TABLE1

        out.print(totalReads);out.print('\t');  // total reads
        //todo add pairing info
        out.print(altAlignments);out.print('\t');  
        out.print(failed);out.print('\t');  //
        out.print(readLength);  //
        //out.print(unique);out.print('\t');      // Unique (non duplicate) reads
        //out.print(duplicates);out.print('\t');  // duplicate reads
        //out.print((float)duplicates / (float) totalReads); // rate of duplicates across all reads
        //todo add estimated library size calculation here
        out.println();

        // TABLE 2
        out.println("Mapped\tMapping Rate\tMapped Unique\tMapped Unique Rate of Total\tUnique Rate of Mapped\tDuplication Rate of Mapped\tBase Mismatch Rate");

        out.print(mappedReads); out.print('\t');
        out.print((float)mappedReads/(float)totalReads);out.print('\t'); // map rate is total mapped divided by total reads
        out.print(mappedUniqueReads); out.print('\t');
        out.print((float)mappedUniqueReads/(float)totalReads); out.print('\t');// map rate is total mapped divided by total reads
        out.print((float)mappedUniqueReads/(float)mappedReads); out.print('\t');// map rate is total mapped divided by total reads
        out.print((float)duplicates/(float)mappedReads); out.print('\t');// map rate is total mapped divided by total reads
        out.print((float)mismatchBases/(float)totalBases); // map rate is total mapped divided by total reads
        out.println();

        // TABLE 3 - pairing information
        out.println("Mapped Pairs\tUnpaired Reads\tEnd 1 Mapping Rate\tEnd 2 Mapping Rate\tEnd 1 Mismatch Rate\tEnd 2 Mismatch Rate\tFragment Length Mean\tFragment Length StdDev\tChimeric Pairs");
        out.print(totalMappedPairs); out.print('\t');
        out.print(unpairedReads); out.print('\t');
        out.print((float)end1Mapped/((float)totalReads/2f)); out.print('\t');
        out.print((float)end2Mapped/((float)totalReads/2f)); out.print('\t');
        out.print((float)end1MismatchBases/(float)end1TotalBases);  out.print('\t');// map rate is total mapped divided by total reads
        out.print((float)end2MismatchBases/(float)end2TotalBases);  out.print('\t');// map rate is total mapped divided by total reads
        float mean = getFragmentLengthMean();
        out.print((int)mean); out.print('\t');
        //out.print(getFragmentLengthMedian()); out.print('\t');
        out.print((int)getFragmentLengthStdDev(mean)); out.print('\t');
        out.print(chimericPairs);
        out.println();
        
        out.println("Intragenic Rate\tExonic Rate\tIntronic Rate\tIntergenic Rate\tSplit Reads\tExpression Profiling Efficiency\tTranscripts Detected\tGenes Detected"); // TABLE 4
        out.print((float)intragenic/(float)mappedReads);out.print('\t');
        out.print((float)exonic/(float)mappedReads);out.print('\t');
        out.print((float)intronOrUTR/(float)mappedReads);out.print('\t');
        out.print((float)intergenic/(float)mappedReads);out.print('\t');
        out.print(splitReads);out.print('\t');
        out.print((float)exonic/(float)totalReads);out.print('\t');
        out.print(transcriptLog.size());out.print('\t');
        out.print(geneLog.size());
        out.println();

        out.println("End 1 Sense\tEnd 1 Antisense\tEnd 2 Sense\tEnd 2 Antisense\tEnd 1 % Sense\tEnd 2 % Sense"); // table 5
        out.print(end1Sense);out.print('\t');
        out.print(end1Antisense);out.print('\t');
        out.print(end2Sense);out.print('\t');
        out.print(end2Antisense);out.print('\t');
        out.print(100f*(float)end1Sense/(float)(end1Antisense+end1Sense));out.print("\t");
        out.print(100f*(float)end2Sense/(float)(end2Sense+end2Antisense));

        out.println();

        try {
            chOut.close();
        } catch(IOException e) {
            throw new RuntimeException(e.getMessage());
        }

        if (fileOk && OUT_FILE != null) {
            out.close();
        } else if (!fileOk) {
            throw new RuntimeException("Failed to write to " + OUT_FILE);
        }
    }


    public static final String TABLE1_FORMAT = "0000"; // is decimal?
    public static final String TABLE2_FORMAT = "0101111";
    public static final String TABLE3_FORMAT = "001111000"; // PAIRS
    public static final String TABLE4_FORMAT = "11110100";
    public static final String TABLE5_FORMAT = "000011";

    
    private float getFragmentLengthMean() {
        float sum = 0;
        float n = 0;

        for (Integer length:  this.fragmentLengths.getKeySet()){
            int count = this.fragmentLengths.getCount(length);
            for (int i = 0 ; i < count ; i ++){
                sum += length;
                n++;
            }
        }
        return sum / n;
    }

    private float getFragmentLengthStdDev(float mean) {
        float sum = 0;
        float n = 0;

        for (Integer length:  this.fragmentLengths.getKeySet()){
            int count = this.fragmentLengths.getCount(length);
            for (int i = 0 ; i < count ; i ++){
                float diff = length - mean;
                sum += (diff * diff);
                n++;
            }
        }
        return (float)Math.sqrt((double)sum / (n-1f));  
    }

}

