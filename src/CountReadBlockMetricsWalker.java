package org.broadinstitute.cga.rnaseq.gatk;

import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;

import java.util.ArrayList;


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
@Requires({DataSource.READS, DataSource.REFERENCE_BASES})
public class CountReadBlockMetricsWalker extends CountReadMetricsWalker {
//    @Argument(fullName = "outfile_metrics", shortName = "OM", doc="The destination file for the metrics", required=true)
//    private String OUT_FILE = null;
//
//    @Argument(fullName="refseq", shortName="refseq",
//            doc="Name of RefSeq transcript annotation file. ", required=false)
//    String RefseqFileName = null;



    @Override
    protected void makeRefSeqDerviedCounts(SAMRecord read, ArrayList<RefSeqFeature> refSeqs, GenomeLoc readLoc) {
        //Performance forRefs = new Performance("RefSeq get location loop",Performance.Resolution.milliseconds);
        boolean wasIntragenic = false;
        boolean wasExonic = false;
        boolean wasIntron = false;
        //boolean wasCoding = false;

        ArrayList<GenomeLoc> blocks = getBlocksForRead(read);
        for (RefSeqFeature refSeq: refSeqs){
            for (GenomeLoc blockLoc : blocks){
        //                        System.out.println("\trefseq:"+refSeq.getGeneName());
        //                        System.out.println("\trefseq:" + refSeq.getTranscriptId());
        //                        System.out.println("\trefseq:"+refSeqLoc);


                //GenomeLoc readLoc = ref.getLocus();
                // sanity check:
                if (!refSeq.overlapsP(blockLoc)){
                    System.out.println("WHY DOESN'T THIS READ OVERLAP IT'S OWN ROD ENTRY !?!??!!?");
                    System.out.println("\tReadLoc:"+readLoc);
                    System.out.println("\trefSeqLoc:"+refSeq.getLocation());
                    System.out.println("\tblockSeqLoc:"+blockLoc);
                    System.out.println("\tRead: " + read);
                    System.out.println("\tTransc:" + refSeq.getTranscriptId());
                }
                //intragenic++;// we are in a gene, otherwise, we would be in the ROD
                wasIntragenic = true;
        //                        System.out.println("\tIs Intragenic ...");
                //Performance overlapps = new Performance("Overlap perf: " , Performance.Resolution.milliseconds);



                if (refSeq.overlapsExonP(blockLoc)){ // returns true if any of the exons overlap this location
        //                            System.out.println("\tIs Exonic");
                    //exonic++;
                    wasExonic=true;

                }else{
                    // if we're in a gene but not in an exon ... it must be an intron? it could be UTR :(
                    //intronOrUTR++;
                    wasIntron = true;
        //                            System.out.println("\tNOT Exonic");
                }

//                if (refSeq.overlapsCodingP(blockLoc)){
//                    //coding++;
//                    //wasCoding = true;
//        //                            System.out.println("\tIs Coding");
//                }

                //overlapps.outputIfTookTime();
            }// end of blocks

        }// end of ROD entries
        if (wasIntragenic) intragenic++;
        if (wasExonic) exonic++;
        if (wasIntron) intronOrUTR++;
        //if (wasCoding) coding++;




    }

    /**
     *  Breaks the read into intervals to exclude large inserts (typically introns).
     * @param read
     * @return a list of intervals, each of which corresponds to a segment of the read
     */
    protected ArrayList<GenomeLoc> getBlocksForRead(SAMRecord read) {
        ArrayList<GenomeLoc> blocks = new ArrayList<GenomeLoc>();

        int start  = read.getAlignmentStart();
        String contig = read.getReferenceName();


        for (final CigarElement e : read.getCigar().getCigarElements()) {
            switch (e.getOperator()) {
                case H : break; // ignore hard clips
                case P : break; // ignore pads
                case S : break; // soft clip read bases
                case N : start += e.getLength(); break;  // advancing start position beyond unaligned segments
                case D : start += e.getLength(); break;
                case I : break;

                case M :
                case EQ :
                case X :
                    // here we are in an aligned block and can create a Locus
                    int length = e.getLength();
                    GenomeLoc block  = getToolkit().getGenomeLocParser().createGenomeLoc(contig, start, start + length -1);
                    blocks.add(block);
                    start  += length;
                    //todo keep track of total length?
                    break;
                default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
            }
        }
        return blocks;
    }

    @Override
    public void onTraversalDone(Integer result) {
        super.onTraversalDone(result); // creates a tab delmited file with the count results

    }

  
}


