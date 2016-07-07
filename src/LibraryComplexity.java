package org.broadinstitute.cga.picardbased;

import net.sf.picard.sam.SortSam;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.cga.tools.IOTools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Created by the Cancer Genome Analysis Group at the Broad Institute.
 * Author: David S. DeLuca
 * User: ddeluca
 * Date: 7/11/11
 * Time: 9:54 AM
 * Developed as part of the RNA-seq analysis efforts of the Broad Institute
 */
public class LibraryComplexity {

    LibraryComplexityResult result = new LibraryComplexityResult();
    boolean isPaired = true;

    public LibraryComplexity(boolean isPaired) {
        this.isPaired = isPaired;
    }

    public LibraryComplexity(LibraryComplexityResult result) {
        result.clearLC(); // the idea is that we can add on to the old result in case we want to merge bams
    }


    public LibraryComplexityResult countFragmentsV2(String inputFile) {
        if( !this.isPaired) return this.countFragmentsSingle(inputFile);

        final SAMFileReader referenceSam = new SAMFileReader(new File(inputFile));
        referenceSam.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);

        for (final SAMRecord read : referenceSam) {
            if (!read.getNotPrimaryAlignmentFlag()) {
                if(!read.getReadPairedFlag()) {
                    result.unPairedReads++;
                } else {
                    //paired read. only count end1
                    if (read.getReadPairedFlag() && read.getFirstOfPairFlag()) {
                        if (read.getReadFailsVendorQualityCheckFlag()) {  // failed vendor QC
                            result.vendorFails++;
                        } else if (read.getReadUnmappedFlag() || read.getMateUnmappedFlag()) {  // pair not aligned
                            result.unalignedPairs++;
                        } else {  // valid aligned pair
                            // check duplication status:
                            boolean unique = !read.getDuplicateReadFlag();
                            if (unique) {
                                result.uniqueFragments++;
                            } else { // duplicate is true
                                result.duplicateFragements++;
                            }
                        }
                    }
                }
            }
        }
        referenceSam.close();
        return result;
    }


    /**
     * Use this method if the bam needs to be sorted by name before calculating library size
     *
     * @param inputFile     bam
     * @param outputFile  can be null if reads are not paired, otherwise this should be the location where a
     *                      temporary sorted bam can be created
     * @return
     */
    public LibraryComplexityResult countFragments(String inputFile, String outputFile) {

        if (isPaired) {
            String [] sortArgs = {"INPUT="+inputFile, "OUTPUT="+outputFile, "SO=queryname", "VALIDATION_STRINGENCY=LENIENT"};
            SortSam sort = new SortSam();
            sort.instanceMain(sortArgs);
        } else {  // no need to sort if not paired
            outputFile = inputFile;
        }
        System.out.println("Counting fragments for library size estimation");
        return computeComplexity(outputFile);
    }


    public LibraryComplexityResult computeComplexity(String inputFile) {
        if (isPaired) return countFragmentsPaired(inputFile);
        else return countFragmentsSingle(inputFile);
    }


    /**
     * use this method if the bam is already sorted by name
     * @param inputFile
     * @return
     */
    private LibraryComplexityResult countFragmentsPaired(String inputFile) {
        final SAMFileReader referenceSam = new SAMFileReader(new File(inputFile));
        referenceSam.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
        SAMRecord r1 = null;
        int count = 0;
        for (final SAMRecord read : referenceSam) {
            count++;
            if (!read.getNotPrimaryAlignmentFlag()) {
                if (r1 == null) {
                    r1 = read;
                } else {
                    // could be second of pair
                    if (read.getReadName().equals(r1.getReadName())) {
                        // alignment issues:
                        if (r1.getReadFailsVendorQualityCheckFlag() || read.getReadFailsVendorQualityCheckFlag()) {
                            result.vendorFails++;
                        } else if (r1.getReadUnmappedFlag() || read.getReadUnmappedFlag()) {
                            result.unalignedPairs++;
                        } else {
                            boolean unique = !r1.getDuplicateReadFlag() && !read.getDuplicateReadFlag();
                            boolean duplicate = r1.getDuplicateReadFlag() && read.getDuplicateReadFlag();
                            if (!unique && !duplicate) {
                                result.inconsistency++;
                            } else {
                                if (unique) {
                                    result.uniqueFragments++;
                                } else { // duplicate is true
                                    result.duplicateFragements++;
                                }
                            }
                        }
                        r1 = null;  // r1 was evaluated and is now cleared
                    } else {
                        result.unPairedReads++;
                        r1 = read; // read didn't match r1, but maybe the next read will match this read
                    }
                }
            }
        }
        referenceSam.close();
        return result;
    }


    /**
     * use this metho d if the bam is already sorted by name
     * @param inputFile
     * @return
     */
    private LibraryComplexityResult countFragmentsSingle(String inputFile) {
        final SAMFileReader referenceSam = new SAMFileReader(new File(inputFile));
        referenceSam.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);

        for (final SAMRecord read : referenceSam) {
            if (!read.getNotPrimaryAlignmentFlag()) {
                result.unPairedReads++;
                if (read.getReadFailsVendorQualityCheckFlag()) {
                    result.vendorFails++;
                } else if (read.getReadUnmappedFlag()) {
                    result.unalignedPairs++;
                } else {
                    if (!read.getDuplicateReadFlag()) {
                        result.uniqueFragments++;
                    } else {  // duplicate is true
                        result.duplicateFragements++;
                    }
                }
            }
        }
        referenceSam.close();
        return result;
    }


    public void toFile(String filename) throws IOException {
        BufferedWriter out = new BufferedWriter(new FileWriter(filename));
        out.write(result.toString());
        out.close();
    }


    public static String getFromFile (String filename)throws IOException {
        try {
            return IOTools.fileToListofArrays(filename,"\\t").get(4)[1];
        } catch (RuntimeException e) {
            System.err.println("Error parsing this file: " + filename);
            System.err.println("File content:");
            System.err.println(IOTools.fileToString(filename));
            throw e;
        }
    }


    public static class LibraryComplexityResult {
        int uniqueFragments = 0;
        int duplicateFragements = 0;
        int unalignedPairs = 0;
        int vendorFails = 0;
        int estimatedLibrarySize = -1;
        int inconsistency =0;
        public int unPairedReads = 0;

        private int estimateLibrarySize() {
            if (duplicateFragements == 0) return 0;
            double n = uniqueFragments+ duplicateFragements;

            int minX = 1;
            int minErr = Integer.MAX_VALUE;

            for (double x = uniqueFragments; x < 1e9; x++) {
                int u = lander(n,x);
                int err = Math.abs(u-uniqueFragments);
                if (err < minErr) {
                    minErr = err;
                    minX = (int)x;
                }
            }
            return minX;
        }


        /**
         * returns unique number of reads
         *
         * @param n
         * @param x
         * @return
         */
        private static int lander(double n, double x) {
            return (int) (
                    (x) * (1d - Math.exp(-1d * n / x))
            );
        }


        public int getEstimatedLibrarySize() {
            if(estimatedLibrarySize ==-1) {
                estimatedLibrarySize = estimateLibrarySize();
            }
            return estimatedLibrarySize;
        }


        public void clearLC() {
            this.estimatedLibrarySize = -1;
        }


        public String toString() {
            return "Unique\t"+uniqueFragments +"\n" +
                    "Duplicate\t"+duplicateFragements+"\n" +
                    "Unaligned\t"+unalignedPairs+"\n" +
                    "Vendor Failed\t"+vendorFails+"\n"+
                    "Estimated Library Size\t"+getEstimatedLibrarySize()+"\n"+
                    "Duplicate Inconsistency\t"+inconsistency +"\n"+
                    "Unpaired Reads\t"+unPairedReads+"\n";
        }
    }

}
