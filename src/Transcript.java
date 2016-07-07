package org.broadinstitute.cga.rnaseq;

import org.broadinstitute.cga.tools.DynamicStreams;
import org.broadinstitute.cga.tools.ObjectCounter;
import org.broadinstitute.cga.tools.Performance;

import java.io.*;
import java.util.*;

public class Transcript {

    // public static void main(String[] args) {
    //     try {
    //         TranscriptList transcripts = TranscriptList.loadGTF("/Users/ddeluca/Documents/RNA-seq/metrics/test/chr18.gtf",null);
    //         for (Transcript t: transcripts) {
    //             if (t.getTranscriptId().equals("NM_001071")) {
    //                 System.out.println(t.toString());
    //
    //                 if (t.getStart() >=0){
    //                     System.out.println("It should be outputed");
    //                 }else{
    //                     System.out.println("it wasn't outputted");
    //                 }
    //             }
    //         }
    //         System.out.println("Done");
    //     }catch(Exception e){
    //         e.printStackTrace();
    //     }
    // }
    
    
//	static final int NAME_COL = 0;     // index values for the expression.txt file
//	static final int ID_COL = 1;
//	static final int BIOTYPE_COL = 2;
//	static final int LENGTH_COL = 3;
//	static final int TOTCOV_COL = 4;
//	static final int AV_COL = 5;
//	static final int STD_COL = 6;

    String transId;
    //String[] expressionData; //Gene_Symbol     Ensembl_Transcript_ID   Ensembl_biotype Transcript_Size Total_Transcript_Coverage       Average_Transcript_Coverage     Transcript_Coverage_StdDev

    String chr = null;

    ArrayList<int[]> exons = new ArrayList<int[]>(); // int[0] is start, int[1] is stop positions
    private float expressionLevel;
    private Annotations annotation = null;
    private DoCMetrics docMetrics;
    private int length = -1;

    private int transcriptStart =-1;
    private int transcriptStop =-1;

    private int cdsStart =-1;
    private int cdsStop =-1;

    Transcript (String id, float expressionLevel) {
        this.transId = id;
        this.expressionLevel = expressionLevel;
    }


    /**
     * returns the cDNA positions of the exon boundries
     * @return
     */
    public Collection<Integer> getLocalExonBoundries() {
        ArrayList<Integer> newExons = new ArrayList<Integer>();

        int lastPos = 0; // the position of the end of the exon
        // The boundaries are the end of the exons, so we can skip the last exon
        for (int i = 0; i < (exons.size()-1); i++){
            int[] ex = exons.get(i);
            int exonSize = ex[1] - ex[0]+1;
            lastPos += exonSize;
            newExons.add(lastPos);
        }
        return newExons;
    }


    public ArrayList<int[]> getExons() {
        return exons;
    }


    public String getTranscriptId() {
        return this.transId;
    }


    public boolean hasPosition(String chr2, int pos) {
        if (!this.chr.equals(chr2)) {
            return false;
        }
        for (int[] exon: exons) {
            if (pos >= exon[0] && pos <=exon[1]) {
                return true;
            }
        }
        return false;
    }


    public float getExpression() {
        return this.expressionLevel;
    }


    void setExonList(String list) {
        this.chr = list.substring(0,list.indexOf(':'));

        this.exons = new ArrayList<int[]>();
        StringTokenizer toks = new StringTokenizer(list,",");
        while (toks.hasMoreElements()){
            String intervalStr = toks.nextToken();
            this.exons.add(parseInterval(intervalStr));
        }
    }


    void setCoordinates(String chrom,  ArrayList<int[]> exons_) {
        this.chr = chrom;
        this.exons = exons_;
    }


    public static int [] parseInterval(String intervalStr) {
        int [] interval = new int[2];
        interval[0] = Integer.parseInt(intervalStr.substring(intervalStr.indexOf(':')+1, intervalStr.indexOf('-')));
        interval[1] = Integer.parseInt(intervalStr.substring(intervalStr.indexOf('-')+1));
        return interval;
    }


    public String getChr() {
        return this.chr;
    }


    public void setAnnotation(Annotations annotation) {
        this.annotation = annotation;
    }


    public Annotations getAnnotation() {
        return annotation;
    }


    public void setDocMetrics(DoCMetrics docMetrics) {
        this.docMetrics = docMetrics;
    }


    public DoCMetrics getDoCMetrics() {
        return docMetrics;
    }


    @Override
    public String toString() {
        String exons = this.getExonList();
        return this.transId + " " + this.chr+":"+exons + " " + annotation.toString();
    }


    /**
     *
     * @param filename
     * @return
     * @throws IOException
     */
    public static Set<String> loadGTFIds(String filename) throws IOException {

        HashSet<String> ids = new HashSet<String>();
        BufferedReader in = DynamicStreams.getDynamicReader(filename);

        String line = in.readLine(); // no header
        while (line != null){
            if (!line.startsWith("#")){
                try {
                    String [] split = line.split("\\t");
                    String[] attributeList = split[8].split(";"); //todo make this more intelligent so order of attribute list doesn't matter
                    String thisId = attributeList[1].replace("transcript_id \"","").replace("\"","").trim(); // extract the id
                    ids.add(thisId);
                } catch (RuntimeException e) {
                    System.out.println("Error parsing line:\t" + line);
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    throw e;
                }
            }
            line = in.readLine();
        }
        in.close();
        return ids;
    }


    public void addExon(int startPos, int stopPos) {
//        int exonIndex = exonNum -1;
//        if (exons == null){
//            exons = new ArrayList<int[]>();
//        }
//
//        int[] e = new int[2];
//        e[0] = startPos;
//        e[1] = stopPos;
//
//        while (exonIndex >= exons.size()){
//            exons.add(null);
//        }
//        exons.set(exonIndex,e);
        // don't add exons without length:

        if ((stopPos-startPos)<=0) return;

        int[] e = new int[2];
        e[0] = startPos; e[1] = stopPos;

        if (exons == null) {
            exons = new ArrayList<int[]>();
            exons.add(e);
            return;
        }

        int index = 0;
        for (; index < exons.size(); index++) {
            if (exons.get(index)[0] > startPos) {
                // we've gone too far.
                break;
            }
        }

        try {
            exons.add(index,e);
        } catch (RuntimeException exc) {
            System.out.println("Index was " + index);
            System.out.println("Length was " + exons.size());
            exc.printStackTrace();
            throw exc;
        }
    }


    /**
     * Produces a string of all the exons as comma separated intervals
     * @return
     */
    public String getExonList() {
        String list = "";

        if (exons == null || exons.size() < 1) {
            return "";
        }

        int[] first = this.exons.get(0);
        list+=first[0]+"-"+first[1];

        for (int i = 1 ; i < exons.size(); i++) {
            int[] e = this.exons.get(i);
            list+=","+e[0]+"-"+e[1];
        }
        return list;
    }


    public static void createIntervalFile(ArrayList<Transcript> transcripts, String intervalFile) throws IOException{
        Performance perf = new Performance("Transcript objects to interval list conversion");
        System.out.println("Writing intervals from transcript objects");

        BufferedWriter out = new BufferedWriter(new FileWriter(intervalFile));

        for (Transcript t:transcripts){
            ArrayList<int[]>  exons = t.getExons();
            if (exons != null){
                out.write(t.getIntervalList());
            }
        }
        out.close();
        System.out.println(perf);
    }


    public String getIntervalList() {
        StringBuilder str = new StringBuilder();
        for (int[] exon: exons){
            str.append(this.getChr()).append(":").append(exon[0]).append("-").append(exon[1]).append("\n");
        }
        return str.toString();
    }


    /**
     * Creates an ordered refgene file from the list of transcript objects. The order is sorted first by the given contig order
     * and then by the start position of the transcript
     *
     * Also modifies the chromosome contig names in the transcript objects to contain those of the reference
     *
     * FOR SOME STUPID REASON UCSC'S REFGEN FORMAT USES 0-BASED START POSITIONS BUT 1-BASED END POSITIONS!!!!!!!
     *
     * @param orderedContigNames
     * @param transcripts
     * @param refGen
     * @throws IOException
     */
    public static void transcriptsToRefGeneFile(ArrayList<String> orderedContigNames,ArrayList<Transcript> transcripts, String refGen) throws IOException{
        Performance perf = new Performance("Transcript objects to RefGen format");
        BufferedWriter out = new BufferedWriter(new FileWriter(refGen));

        for (String contig: orderedContigNames) {
            // for a contig, we must first filter, sort and then output

            // first sort for the transcripts in this contig
            TreeSet<Transcript> orderedTranscripts = new TreeSet<Transcript>(new TranscriptComparator());
            boolean contigHasTranscripts = false;
            String contigAlt;
            if (!contig.contains("chr")) {
                contigAlt ="chr"+contig;
            } else {
                contigAlt = contig.replace("chr","");
            }

            for (Transcript t:transcripts) {
                if (t.getChr().equals(contig) || t.getChr().equals(contigAlt)) {
                    t.setChr(contig);
                    orderedTranscripts.add(t);
                    contigHasTranscripts = true;
                }
            }

            // now output the transcripts in order for this contig
            for (Transcript t:orderedTranscripts) {
                if (t.getStart() >=0) { // we don't want to be outputting transcripts with no start
                    out.write("0"); out.write('\t'); // bins ... hopefully they don't matter
                    out.write(t.getTranscriptId()); out.write('\t'); // name
                    out.write(t.getChr()); out.write('\t'); // chrom
                    out.write(t.getAnnotation().getStrand()); out.write('\t'); // strand

                    out.write(""+(t.getStart()-1)); out.write('\t'); // start
                    out.write(""+t.getStop()); out.write('\t'); // stop

                    out.write(""+(t.getCdsStart()-1)); out.write('\t'); // start
                    out.write(""+t.getCdsStop()); out.write('\t'); // stop

                    out.write(""+t.getExons().size()); out.write('\t'); // exonCount

                    //exon starts:
                    for (int[] exon: t.getExons()){
                        out.write(""+(exon[0]-1)+",");
                    }
                    out.write("\t");

                    //exon stops:
                    for (int[] exon: t.getExons()){
                        out.write(""+exon[1]+",");
                    }

                    // score
                    out.write("\t0\t");

                    if(t.getAnnotation().getName() != null) {
                        out.write(t.getAnnotation().getName()); // gene name
                    }
                    out.write('\t');

                    //todo cdsStartStat and  cdsEndStat
                    out.write("unk\tunk\t"); // should be changed to cmpl if it is actually known

                    //exon frames:
                    for (int i = 0 ; i < t.getExons().size(); i++) {
                        out.write("-1,");
                    }

                    out.write('\n');
                }
            }
        }
        out.close();
        System.out.println(perf);
    }

    private void setChr(String contig) {
        this.chr = contig;
    }

    private int getCdsStop() {
        if (cdsStop < 0) return transcriptStop;
        return cdsStop;
    }

    private int getCdsStart() {
        if (cdsStart < 0) return transcriptStart;
        return cdsStart;
    }

    public void setStop(int stop) {
        this.transcriptStop = stop;
    }

    public void setStart(int start) {
        this.transcriptStart = start;
    }


    /**
     * This method dopies a GTF and selects only those transcripts with ids found in the selectedIds Set
     *
     *
     * @param sourceGTF
     * @param destFile
     * @param selectedIds
     * @throws IOException
     */
    public static void copyFiltered(String sourceGTF, String destFile, Collection<String> selectedIds) throws IOException {

        BufferedWriter out = new BufferedWriter(new FileWriter(destFile));
        BufferedReader in = DynamicStreams.getDynamicReader(sourceGTF);

        String line = in.readLine(); // no header
        while (line != null){
            if (!line.startsWith("#")) {
                try {
                    String [] split = line.split("\\t");
                    String[] attributeList = split[8].split(";"); //todo make this more intelligent so order of attribute list doesn't matter

                    String thisId = attributeList[1].replace("transcript_id \"","").replace("\"","").trim(); // extract the id

                    if (selectedIds.contains(thisId)){
                        out.write(line); out.write('\n');
                    }
                } catch (RuntimeException e) {
                    System.out.println("Error parsing line:\t" + line);
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    throw e;
                }
            }
            line = in.readLine();
        }
        in.close();
        out.close();
    }

    public void setGC(Float val) {
        this.annotation.setGCContent(val);
    }

    public float getGCContent(){
        return this.annotation.getGCContent();
    }

    public ArrayList<Integer> getCoverage() {
        return this.docMetrics.getCoverage();
    }

    public void merge(Transcript t) {
        ArrayList<int[]> newExons = new ArrayList<int[]>();
        int l = t.exons.size() + this.exons.size();
        int[] starts = new int[t.exons.size() + this.exons.size()];
        int i = 0;
        
        for (int j = 0 ; j < t.exons.size(); j++){
            starts[i] = t.exons.get(j)[0];
            i++;
        }
        for (int j = 0 ; j < this.exons.size(); j++){
            starts[i] = this.exons.get(j)[0];
            i++;
        }
        
        int[] stops = new int[t.exons.size() + this.exons.size()];
        i=0;
        for (int j = 0 ; j < t.exons.size(); j++){
            stops[i] = t.exons.get(j)[1];
            i++;
        }
        for (int j = 0 ; j < this.exons.size(); j++){
            starts[i] = this.exons.get(j)[1];
            i++;
        }
        
        Arrays.sort(starts);
        Arrays.sort(stops);

        int startI =0;
        int stopI = 0;
        int thisPos = starts[startI];
        startI++;
        do {
            // if we have another start before the next stop
            if (starts[startI] < stops[stopI]){
                startI++;
            } else {
                // the next start is after this next stop, so we conclude the exon
                int[] newExon = new int[2];
                newExon[0] = thisPos;
                newExon[1] = stops[stopI];
                stopI++;
                newExons.add(newExon);
                thisPos = starts[startI];
                startI++;
            }
        } while (startI < starts.length && stopI < stops.length);
    }


    public static class TranscriptComparator implements Comparator<Transcript> {

        @Override
        public int compare(Transcript t1, Transcript t2) {
            if (!t1.getChr().equals(t2.getChr())) {
                return t1.getChr().compareTo(t2.getChr());
            }

            Integer i1 = t1.getStart();
            Integer i2 = t2.getStart();

            int result = i1.compareTo(i2);
            if (result == 0){
                return t1.getTranscriptId().compareTo(t2.getTranscriptId());
            }
            return result;
        }
    }

    private int getStart() {
        // -1 means uninitialized
        if (transcriptStart == -1){
            if (exons.size() > 0) {
                int [] e = exons.get(0);
                transcriptStart = e[0];
            } else {
                // no hope
                transcriptStart = -2;
            }

        }
        return transcriptStart;
    }

    private int getStop() {
        if (transcriptStop == -1){
            if (exons.size() > 0) {
                int [] e = exons.get(exons.size()-1);
                transcriptStop = e[1];
            } else {
                // no hope
                transcriptStop = -2;
            }
        }
        return transcriptStop;
    }

    /**
     * returns the length of the transcript based on all the exon lengths
     * @return
     */
    public int getLength() {
        if (length < 0){ // calculate length from exons on first call
            length = 0;
            if (this.exons != null) {
                for (int[] e: exons) {
                    length += e[1] - e[0] +1;
                }
            }
        }
        return length;  //To change body of created methods use File | Settings | File Templates.
    }


    /**
     * RPKM: fragments per kilobase of exon per million fragments mapped
     *
     * @param exonReads
     * @param exonLength
     * @param mappedReads
     */
    public static float getRPKM(int exonReads, int exonLength, int mappedReads) {
        float kbExon = ((float)exonLength) / 1000f;
        float milFragMapped = ((float)mappedReads) / 1E6f;
        return ((float)exonReads) / kbExon / milFragMapped;
    }


    /**
     *
     * @param infile
     * @return
     * @throws IOException
     */
    public static TranscriptList readTranscriptsFromList(String infile) throws IOException{
        BufferedReader in = new BufferedReader(new FileReader(infile));
        TranscriptList trans = new TranscriptList();

        String line = in.readLine();
        while (line!=null){
            trans.add(new Transcript(line.trim(), 0));
            line = in.readLine();
        }
        in.close();
        return trans;
    }


    public int binGaps(ObjectCounter gapLengthCounter, int MAX_BIN){
        int thisGap =-1;
        ArrayList<Integer> cov = this.getCoverage();
        int longGaps =0;
        for (Integer c: cov){
            if (c <= TranscriptDoCResult.MAX_COV_FOR_GAP){
                // we are in a gap
                if (thisGap == -1) {
                    // we are starting a new gap
                    thisGap =0;
                }
                thisGap++; // extend the gap
            } else {
                // no gap.
                // if we are ending a gap:
                if (thisGap > 0) {
                    // save last gap
                    if (thisGap >= TranscriptDoCResult.MIN_GAP_SIZE) {
                        if (thisGap >= MAX_BIN) {
                            longGaps++;
                        } else {
                            gapLengthCounter.add(thisGap);
                        }
                    }
                    thisGap = -1;
                }
            }
        }
        // check whether the end is a gap
        if (thisGap > 0) {
            if (thisGap >= TranscriptDoCResult.MIN_GAP_SIZE) {
                if (thisGap >= MAX_BIN) {
                    longGaps++;
                } else {
                    gapLengthCounter.add(thisGap);
                }
            }
        }
        return longGaps;
    }

}
