package org.broadinstitute.cga.rnaseq.gatk;

import com.sun.corba.se.impl.orbutil.ObjectUtility;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.DailyRollingFileAppender;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.RODMetaDataContainer;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.SeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.awt.image.ImagingOpException;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;


/**
 * Created by the Broad Institute, Cancer Genome Analysis Group
 * Author: David S. DeLuca
 * User: ddeluca
 * Date: 5/19/11
 * Time: 10:30 AM
 *
 * This adaptation of the CountReadWalker that outputs the total count into a file at the end
 * (total reads <t a b > unique reads
 *
 */
@Requires({DataSource.READS, DataSource.REFERENCE_BASES})
public class OutputCountReadsWalker extends ReadWalker<Integer,Integer> {
    @Argument(fullName = "outfile_readCounts", shortName = "OC", doc="The destination file for the counts", required=true)
    protected String OUT_FILE = null;

    int totalReads =0;
    int uniqueReads =0;


    @Override
    public Integer map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker) {

        
        if (!read.getNotPrimaryAlignmentFlag() && !read.getReadFailsVendorQualityCheckFlag() ){
            totalReads++;
            if (!read.getDuplicateReadFlag()){
                uniqueReads++;
            }
        }
        return 1;

    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    @Override
    public void onTraversalDone(Integer result) {
        super.onTraversalDone(result);

        PrintWriter out;
        boolean fileOk = true; // we're going to do some wierd handling to output to file if it's there. if
        // no file is given, we'll go to stdout or if the file fails to open
        if (OUT_FILE == null){
            out  = new PrintWriter(System.out);
        }else{
            try {
                out = new PrintWriter (new FileWriter(OUT_FILE ));
            } catch (IOException e){
                e.printStackTrace();
                System.out.println("IOException for file: " + OUT_FILE);
                fileOk=false;
                out = new PrintWriter(System.out);
            }
        }

        out.write(""+totalReads+"\t"+uniqueReads);


        if (fileOk && OUT_FILE != null){
            out.close();
        }else if (!fileOk){
            throw new RuntimeException("Failed to write to " + OUT_FILE);
        }
    }


}


