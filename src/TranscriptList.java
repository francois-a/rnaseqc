package org.broadinstitute.cga.rnaseq;

import org.broadinstitute.cga.tools.DynamicStreams;
import org.broadinstitute.cga.tools.Performance;

import java.io.*;
import java.util.*;

/**
 * Created by the Cancer Genome Analysis Group at the Broad Institute.
 * Author: David S. DeLuca
 * User: ddeluca
 * Date: 8/4/11
 * Time: 3:40 PM
 * Developed as part of the RNA-seq analysis efforts of the Broad Institute
 */
public class TranscriptList extends ArrayList<Transcript> {

    public static void main (String[] args) {
        try {
            Performance perf = new Performance("Total Time");

            ArrayList<int[]> exs1 = new ArrayList<int[]>();
            exs1.add(new int[] {1, 5});
            exs1.add(new int[] {10, 15});
            ArrayList<int[]> exs2 = new ArrayList<int[]>();
            exs2.add(new int[] {1, 5});
            exs2.add(new int[] {10, 15});
            exs2.add(new int[] {20, 25});

            ArrayList<int[]> newExons = new ArrayList<int[]>();
            int index1 = 0;
            int index2 = 0;
            int minStart = Math.min(exs1.get(index1)[0],exs2.get(index2)[0]);
            int minStop = Math.min(exs1.get(index1)[1],exs2.get(index2)[1]);

            System.out.println("***");
            for (int[] e: newExons) {
                System.out.println("Ex:\t" + e[0] + "\t" +e[1]);
            }
            System.out.println("Done. " + perf);
        } catch(Exception e) {
            System.out.println("Fatal Error");
            e.printStackTrace();
        }
    }
    
    public static TranscriptList collapse(TranscriptList transcripts) {
        HashMap<String,Transcript> byGeneId = new HashMap<String,Transcript>();
        
        for (Transcript t: transcripts) {
            Transcript gene = byGeneId.get(t.getAnnotation().getGeneId());
            if (gene == null) {
                byGeneId.put(t.getAnnotation().getGeneId(),t);
            } else {
                gene.merge(t);
            }
        }
        System.out.println("Gene Ids\t" + byGeneId.size());
        return null;
    }


    /**
     *
     * @param filename
     * @param transcriptTypeColumn this should have the value, 2, for Ensembl, because Ensembl packs the transcript type into the second column
     * @return
     * @throws IOException
     */
    public static TranscriptList loadGTF(String filename, String transcriptTypeColumn) throws IOException {

        // support ENSEBML's non -spec conforming GTF format (they put transcript_type in the source field)
        TranscriptList  transcripts = new TranscriptList();
        HashMap<String,Transcript> indexed = new HashMap<String,Transcript>();
        BufferedReader in = DynamicStreams.getDynamicReader(filename);

        String line = in.readLine(); // no header
        Transcript currentTranscript = null;

        while (line != null) {
            if (!line.startsWith("#")) {
                try {
                    String [] split = line.split("\\t");
                    //check if the gtf file does not contain any attributes
                    if (split.length < 9) {
                        System.err.println("No attributes were found in the GTF file on line " + line);
                        System.exit(1);
                    }

                    String[] attributeList = split[8].split(";"); //todo make this more intelligent so order of attribute list doesn't matter

                    int index1=-1;
                    for(int a =0; a < attributeList.length;a++) {
                        if(attributeList[a].contains("transcript_id")) {
                            index1 = a;
                        }
                    }
                    if(index1 == -1) {
                        System.err.println("The required transcript_id attribute was not found on line "+ line);
                        System.exit(1);
                    }
                    int transcriptIdIndex = attributeList[index1].indexOf("transcript_id \"");
                    int transcriptIdIndexLength = transcriptIdIndex + 15;

                    String thisId = attributeList[index1].substring(transcriptIdIndexLength,attributeList[index1].length()).trim(); // extract the id
                    thisId = thisId.substring(0, thisId.length()-1);

                    currentTranscript = indexed.get(thisId);
                    if (currentTranscript == null) {
                        currentTranscript = new Transcript(thisId,0);
                        indexed.put(thisId, currentTranscript);
                        transcripts.add(currentTranscript);
                        // add one-time infos:
                        Annotations anno = new Annotations();

                        for (int a =0; a < attributeList.length;a++) {
                            if (attributeList[a].contains("gene_id")) {
                                int index = attributeList[a].indexOf("gene_id \"");
                                int length = index + 9;
                                String attr = attributeList[a].substring(length, attributeList[a].length()).trim(); // extract the id
                                attr = attr.substring(0, attr.length()-1);
                                anno.setGeneId(attr);
                            }
                            if (attributeList[a].contains("gene_name")) {
                                int index = attributeList[a].indexOf("gene_name \"");
                                int length = index + 11;
                                String attr = attributeList[a].substring(length, attributeList[a].length()).trim(); // extract the id
                                attr = attr.substring(0, attr.length()-1);
                                anno.setName(attr); // extract the name
                            }
                            if (attributeList[a].contains("transcript_type")) {
                                int index = attributeList[a].indexOf("transcript_type \"");
                                int length = index + 17;
                                String attr = attributeList[a].substring(length, attributeList[a].length()).trim(); // extract the id
                                attr = attr.substring(0, attr.length()-1);
                                anno.setType(attr); // extract transcript_type
                            }
                        }

                        if (anno.getGeneId() == null) {
                            System.err.println("The required gene_id attribute was not found on line "+ line);
                            System.exit(1);
                        }
                        anno.setStrand(split[6]);

                        // if a column was provided for the trancript type use this
                        // instead of transcript_type attribute
                        if (transcriptTypeColumn!= null) {
                            int column;

                            column = Integer.parseInt(transcriptTypeColumn);
                            //subtract 1 since array indexes are 0 based
                            column = column - 1;
                            //todo allow transcript type column to also be an attribute
                            if (transcriptTypeColumn != null && column > 8) {
                                System.err.println("The transcript type column should be less than or equal to 8");
                                System.exit(1);
                            }
                            anno.setType(split[column]); // extract transcript type from a column in gtf
                        }
                        currentTranscript.setAnnotation(anno);
                        currentTranscript.chr = split[0]; // chromosome
                    } // end of section for new transcripts

                    // now add this line nto transcript
                    if (split[2].equals("exon")) {
                        currentTranscript.addExon( Integer.parseInt(split[3]),Integer.parseInt(split[4]));
                    } else if (split[2].equals("transcript")) {
                        currentTranscript.setStart(Integer.parseInt(split[3]));
                        currentTranscript.setStop(Integer.parseInt(split[4]));
                    }
                } catch (RuntimeException e) {
                    System.out.println("Error parsing line:\t" + line);
                    e.printStackTrace();
                    throw e;
                }
            }
            line = in.readLine();
        }
        in.close();
        return transcripts;
    }


    public void toFile(String outfileName, String endLength, boolean  gc, boolean details) throws IOException {
        BufferedWriter out = new BufferedWriter(new FileWriter(outfileName));

        out.write("Rank\tTranscript ID\tGene Name\tType\tLength\tTotal Cov.\tMean Cov.\tStd Dev\tCV\t5'"+endLength+"Base");
        out.write("\t3' "+endLength+"Base");
        out.write("\t5' "+endLength+"Base Norm");
        out.write("\t3' "+endLength+"Base Norm");
        out.write("\t5':3' Ratio "+endLength+"Base Norm");
        out.write("\tNum. Gaps\tCumul. Gap Length\tGap %");
        if(gc)out.write("\tGC %");
        out.write("\n");

        int rank = 1;
        for (Transcript trans: this) {
            DoCMetrics t = trans.getDoCMetrics();
            out.write(""+ rank );
            rank++;
            if(details) {
                out.write("\t"+ "<a href=" + trans.getTranscriptId() +".html>" + trans.getTranscriptId() +"</a>");
            } else {
                out.write("\t"+ trans.getTranscriptId());
            }
            out.write("\t" +trans.getAnnotation().getName() );	// Gene Name
            out.write("\t" +trans.getAnnotation().getType()); // biotype
            out.write("\t" +trans.getLength() ); // Length
            out.write("\t" +t.getTotalCoverage() ); // Total Cov.
            out.write("\t" +String.format("%.2f",t.getMeanCoverage()) ); // Mean Cov.
            out.write("\t" +String.format("%.2f", t.getStdDev()) ); //Stdev
            out.write("\t" +String.format("%.4f", t.getCV()) ); //CV

            if (endLength.equals("10")) {
                out.write("\t" +t.getFiveEnd10());
            } else if (endLength.equals("50")) {
                out.write("\t" +t.getFiveEnd50());
            } else if (endLength.equals("100")) {
                out.write("\t" +t.getFiveEnd100());
            }
            
            if (endLength.equals("10")) {
                out.write("\t" +t.getThreeEnd10());
            } else if (endLength.equals("50")) {
                out.write("\t" +t.getThreeEnd50());
            } else if (endLength.equals("100")) {
                out.write("\t" +t.getThreeEnd100());
            }

            if (endLength.equals("10")) {
                out.write("\t"+String.format("%.3f", t.getFiveEnd10()/ (t.getMeanCoverage())));
            } else if(endLength.equals("50")) {
                out.write("\t"+String.format("%.3f", t.getFiveEnd50()/(t.getMeanCoverage())));
            } else if (endLength.equals("100")) {
                out.write("\t"+String.format("%.3f", t.getFiveEnd100()/(t.getMeanCoverage())));
            }

            if (endLength.equals("10")) {
                out.write("\t"+ String.format("%.3f", t.getThreeEnd10()/ (t.getMeanCoverage())));
            } else if(endLength.equals("50")) {
                out.write("\t"+String.format("%.3f", t.getThreeEnd50()/(t.getMeanCoverage())));
            } else if (endLength.equals("100")) {
                out.write("\t"+String.format("%.3f", t.getThreeEnd100()/(t.getMeanCoverage())));
            }

            if (endLength.equals("10")) {
                out.write("\t"+String.format("%.3f", (t.getFiveEnd10()/(t.getMeanCoverage())) / (t.getThreeEnd10()/(t.getMeanCoverage()))));
            } else if(endLength.equals("50")) {
                out.write("\t"+String.format("%.3f", (t.getFiveEnd50()/(t.getMeanCoverage())) / (t.getThreeEnd50()/(t.getMeanCoverage()))));
            } else if (endLength.equals("100")) {
                out.write("\t"+String.format("%.3f", (t.getFiveEnd100()/(t.getMeanCoverage())) / (t.getThreeEnd100()/(t.getMeanCoverage()))));
            }

            out.write("\t"+ t.getNumberOfGaps() );
            out.write("\t"+t.getCumulativeGapLength() );
            out.write("\t"+String.format("%.1f",(t.getGapRatio() * 100f)));
            if (gc) {
                out.write("\t"+String.format("%.1f",(trans.getAnnotation().getGCContent() * 100f)));
            }
            out.write("\n");
        }
        out.close();
    }

    /**
     * use map: /Volumes/cga_cga2/berger/RNASeq/refs/Ensembl52.plus.Genome.map
     *
     * each line of this map has an ID[space]comma-delim-chromcoordinates
     *
     * This only creates intervals for those transcripts in our top transcript list
     *
     * @param mapFileName
     * @return
     * @throws IOException
     */
    public Collection<String> loadIntervalsIntoTranscripts(String mapFileName) throws IOException {
        ArrayList<String> intervals= new ArrayList<String>();

        BufferedReader  in = new BufferedReader (new FileReader(mapFileName));

        String line = in.readLine();

        while (line!=null) {
            for (Transcript tran: this) {
                if (line.contains(tran.getTranscriptId())) {
                    intervals.add(line);
                    try {
                        tran.setExonList(line.substring(line.indexOf(' ')+1));
                    } catch (NumberFormatException e) {
                        System.err.println("Could not parse " + line);
                        throw e;
                    }
                }
            }
            line = in.readLine();
        }
        in.close();
        return intervals;
    }


    public void toRRNAIntervalList( String rRNAIntervalFile) throws IOException {
        BufferedWriter out = new BufferedWriter(new FileWriter(rRNAIntervalFile));

        boolean foundrRNA = false;
        for (Transcript t: this) {
            if (t.getAnnotation().getType() != null && t.getAnnotation().getType().contains("rRNA")) {
                if (t.getExons() != null) {
                    out.write(t.getIntervalList());
                    foundrRNA = true;
                }
            }
        }
        out.close();
        if (!foundrRNA) {
            throw new RuntimeException("No rRNA found in GTF transcript_type field");
        }
    }


    public int getTotalCoverage() {
        int tot = 0;
        for(Transcript t: this){
            tot+= t.getDoCMetrics().getTotalCoverage();
        }
        return tot;
    }


    public float getAverageCoverages() {
        int tot = 0; //total length
        for(Transcript t: this){
            tot+= t.getLength();
        }
        return (float)getTotalCoverage() / (float)tot;
    }


    public float getAverageCV() {
        float tot = 0; //total length
        for(Transcript t: this){
            tot+= t.getDoCMetrics().getCV();
        }
        return tot / (float)this.size();
    }


    public float getAvFiveEnd10() {
        float tot = 0;
        for (Transcript t: this) {
            tot+= t.getDoCMetrics().getFiveEnd10Norm();
        }
        return tot / (float)this.size();
    }


    public float getAvFiveEnd50() {
        float tot = 0;
        for (Transcript t: this) {
            tot+= t.getDoCMetrics().getFiveEnd50Norm();
        }
        return tot / (float)this.size();
    }


    public float getAvFiveEnd100() {
        float tot = 0;
        for (Transcript t: this) {
            tot+= t.getDoCMetrics().getFiveEnd100Norm();
        }
        return tot / (float)this.size();
    }


    public float getAvThreeEnd10() {
        float tot = 0;
        for (Transcript t: this) {
            tot+= t.getDoCMetrics().getThreeEnd10Norm();
        }
        return tot / (float)this.size();
    }


    public float getAvThreeEnd50() {
        float tot = 0;
        for (Transcript t: this) {
            tot+= t.getDoCMetrics().getThreeEnd50Norm();
        }
        return tot / (float)this.size();
    }


    public float getAvThreeEnd100() {
        float tot = 0;
        for (Transcript t: this) {
            tot+= t.getDoCMetrics().getThreeEnd100Norm();
        }
        return tot / (float)this.size();
    }


    public int getThreeEndCovCount() {
        int count = 0;
        for (Transcript t: this){
            if (t.getDoCMetrics().has3EndCov()) count++;
        }
        return count;
    }


    public int getFiveEndCovCount() {
        int count = 0;
        for (Transcript t: this) {
            if (t.getDoCMetrics().has5EndCov()) count++;
        }
        return count;
    }


    public int getNumberOfGaps() {
        int tot = 0;
        for (Transcript t: this){
            tot+= t.getDoCMetrics().getNumberOfGaps();
        }
        return tot;
    }


    public int getCumulativeGapLength() {
        int tot = 0;
        for (Transcript t: this){
            tot+= t.getDoCMetrics().getCumulativeGapLength();
        }
        return tot;
    }


    /**
     * Fills the transcripts (i.e. best expressed transcripts), with the chrom and exon coordinate info from the fullTrans list
     * @param fullTrans
     */
    public void loadIntervalsAndAnnotation(ArrayList<Transcript> fullTrans) {
        for (Transcript full: fullTrans) {
            for (Transcript tran: this) {
                if (full.getTranscriptId().equals(tran.getTranscriptId())) {
                    tran.setCoordinates(full.getChr(),full.getExons());
                    tran.setAnnotation(full.getAnnotation());
                }
            }
         }
    }


    public void populateGC(String gcFile) throws IOException {

        //load map
        HashMap<String,Float> gcMap = new HashMap<String,Float>();
        BufferedReader in = new BufferedReader(new FileReader(gcFile));

        String line = in.readLine(); //no header
        while(line !=null) {
            if (!line.trim().equals("")) {
                String [] split = line.split("\\t");
                gcMap.put(split[0],Float.valueOf(split[1]));
            }
            line = in.readLine();
        }
        in.close();

        for(Transcript t: this) {
            t.setGC(gcMap.get(t.getTranscriptId()));
        }
    }


    /**
     * This is the ratio between the cumulative gap length and the cumulative length
     *
     * @return
     */
    public float getGapRatio() {
        float totLength = 0;
        float totalGap = this.getCumulativeGapLength();
        for (Transcript t: this){
            totLength += Integer.valueOf(t.getLength());
        }
        return  totalGap / totLength;
    }


    public void toSummaryFile(String doCTranscriptsSummaryResultsFile) throws IOException {
        BufferedWriter out = new BufferedWriter(new FileWriter(doCTranscriptsSummaryResultsFile));
        out.write(""+this.getAverageCoverages()+"\n");
        out.write(""+this.getAverageCV()+"\n");
        out.write(""+getFiveEndCovCount()+"\n");
        out.write(""+this.getAvFiveEnd10()+"\n");
        out.write(""+this.getAvFiveEnd50()+"\n");
        out.write(""+this.getAvFiveEnd100()+"\n");
        out.write(""+getThreeEndCovCount()+"\n");
        out.write(""+this.getAvThreeEnd10()+"\n");
        out.write(""+this.getAvThreeEnd50()+"\n");
        out.write(""+this.getAvThreeEnd100()+"\n");
        out.write(""+this.getNumberOfGaps()+"\n");
        out.write(""+this.getCumulativeGapLength()+"\n");
        out.write(""+this.getGapRatio()+"\n");
        out.close();
    }
}
