package org.broadinstitute.cga.tools;

import javax.xml.bind.SchemaOutputResolver;
import java.io.*;
import java.util.*;

/**
 * Created by the Cancer Genome Analysis Group at the Broad Institute.
 *
 *
 * Author: David S. DeLuca
 * User: ddeluca
 * Date: 6/22/11
 * Time: 4:45 PM
 * Developed as part of the RNA-seq analysis efforts of the Broad Institute
 */
public class GCTFile {

    String[] ids ;
    String desc[] ;
    String[] samples;
    double [][] data; // [samps,rows]

    int nRows, nSamps;

    public GCTFile(String filename) throws IOException {
        BufferedReader in = new BufferedReader (new FileReader(filename));

        String line = in.readLine(); // that weird #1.2 line
        line = in.readLine(); // line with number of rows and columns
        nRows = Integer.valueOf(line.split("\\t")[0]);
        nSamps = Integer.valueOf(line.split("\\t")[1]);

        ids = new String[nRows];
        desc = new String[nRows];
        samples = new String[nSamps];
        data = new double[nSamps][nRows];

        line = in.readLine(); // headers
        String [] headers = line.split("\\t");
        for (int i = 2; i < headers.length;i++) {
            samples[i-2] = headers[i];
        }

        line = in.readLine(); // actual first data line

        int row = 0;
        while (line!= null){
            String [] split = line.split("\\t"); // form: Accname  desc    sampleName
            ids[row] = split[0];
            desc[row] = split[1];
            for (int i = 2; i < split.length; i++) {
                float val = Float.valueOf(split[i]);
                data[i-2][row]= val;
            }
            line = in.readLine();
            row++;
        }
        in.close();
    }


    public GCTFile(int nRows_, int nSamps_) {
        this.nRows = nRows_;
        this.nSamps = nSamps_;
        ids = new String[nRows];
        desc = new String[nRows];
        samples = new String[nSamps];
        data = new double[nSamps][nRows];
    }


    /**
     * Converts the table into a GCT file. If the table starts with comment (#) rows, these are skipped
     * After the optional comments, a header line is expected
     * The table's first column is the ID and each subsequent column is a sample
     *
     * @param table
     */
    public GCTFile(ArrayList<String[]> table) {
        int start =0;
        // advance passed comments
        while (table.get(start)[0].startsWith("#")) {
            start++;
        }
        String[] header = table.get(start);
        start++;

        nRows = table.size() - start;
        nSamps = header.length -1;

        ids = new String[nRows];
        desc = new String[nRows];
        samples = new String[nSamps];
        data = new double[nSamps][nRows];

        for (int i = 1; i < header.length;i++) {
            samples[i-1] = header[i];
        }

        int row = 0;
        for (int oldRow = start; oldRow < table.size(); oldRow++) {
            String [] split = table.get(oldRow);
            ids[row] = split[0];
            desc[row] = "";
            for (int i = 1; i < split.length; i++) {
                float val = Float.valueOf(split[i]);
                data[i-1][row]= val;
            }
            row++;
        }
    }


    /**
     * Returns a new GCT the combines this object with the parameter object. This object serves as a reference
     * such that the resulting GCT contains the same rows as this object. The values from the parameter objected
     * are quiried using the gene name located in the description field
     *
     *
     * @param gct
     * @return
     * @throws IOException
     */
    public GCTFile combine(GCTFile gct, boolean useDescFromPassedObj, boolean useDescFromThisObj) {
        GCTFile newGCT = new GCTFile(this.nRows , this.nSamps+gct.nSamps);

        int newSampDex = 0;

        for (;newSampDex < this.nSamps; newSampDex++) {
            // add the header
            newGCT.samples[newSampDex]=this.samples[newSampDex];

            // add the data by row
            for (int j = 0; j < nRows; j++) {
                newGCT.data[newSampDex][j] = this.data[newSampDex][j];
                newGCT.ids[j] = this.ids[j];
                newGCT.desc[j] = this.desc[j];
            }
        }

        // now add the info from the second gct
        HashMap<String,double[]> indexed = null;
        if (useDescFromPassedObj) {
            indexed = gct.getIndexedByDesc();
            int i = 0;
            for (;newSampDex < newGCT.nSamps; newSampDex++) {
                // add the header
                newGCT.samples[newSampDex]=gct.samples[i];
                // add the data by row
                for (int j = 0; j < nRows; j++) {
                    double[] row;
                    if(useDescFromThisObj) {
                        row = indexed.get(newGCT.desc[j]);
                    } else {
                        row = indexed.get(newGCT.ids[j]);
                    }
                    if (row==null) {
                        newGCT.data[newSampDex][j] = Float.NaN;
                    } else {
                        newGCT.data[newSampDex][j] = row[i];
                    }
                }
                i++;
            }
        } else {
            indexed = gct.getIndexedById();
            int i = 0;
            for (;newSampDex < newGCT.nSamps; newSampDex++) {
                // add the header
                newGCT.samples[newSampDex]=gct.samples[i];
                // add the data by row
                for (int j = 0; j < nRows; j++) {
                    double[] row;
                    if (useDescFromThisObj) {
                        row = indexed.get(newGCT.desc[j]);
                    } else {
                        row = indexed.get(newGCT.ids[j]);
                    }
                    if (row==null) {
                        newGCT.data[newSampDex][j] = Float.NaN;
                    } else {
                         newGCT.data[newSampDex][j] = row[i];
                    }
                }
                i++;
            }
        }
        return newGCT;
    }


    private HashMap<String, double[]> getIndexedByDesc() {
        HashMap<String,double[]> index = new HashMap<String,double[]>();
        for (int i = 0 ; i< nRows; i++){
            double[] row = new double[nSamps];

            for (int j = 0 ; j < nSamps; j++) {
                row[j]=data[j][i];
            }
            index.put(desc[i],row);
        }
        return index;
    }


    /**
     * the index contains a mapping of the row's id with the row data (a double array)
     * @return
     */
    private HashMap<String, double[]> getIndexedById() {
        HashMap<String,double[]> index = new HashMap<String,double[]>();
        for (int i = 0 ; i< nRows; i++) {
            double[] row = new double[nSamps];
            for (int j = 0 ; j < nSamps; j++) {
                row[j]=data[j][i];
            }
            if (ids[i] == null) {
                throw new RuntimeException ("Warning: When indexing GCT by ID, an ID value of null was found");
            }
            index.put(ids[i],row);
        }
        return index;
    }


    /**
     * This returns the actual array of samples, not a copy, which means that modifications woudl affect hte object
     * @return
     */
    public String[] getSamples() {
        return samples;
    }


    public void toFile(String filename) throws IOException {
        BufferedWriter out = new BufferedWriter(new FileWriter(filename));
        out.write("#1.2\n");

        // no. of rows:
        out.write(""+this.nRows);
        out.write('\t');
        out.write(""+nSamps);
        out.write("\nName\tDescription");
        for (String samp: samples){
            out.write('\t');
            out.write(samp);
        }
        out.write('\n'); // header lines done

        for (int row = 0; row < nRows; row++){
            out.write(ids[row]);
            out.write('\t');
            out.write(desc[row]);

            for (int j = 0 ; j < nSamps; j++){
                out.write('\t');
                out.write(""+data[j][row]);
            }
            out.write('\n');
        }
        out.close();
    }


    /**
     * the finalData objecdt must contain the id and description fields, as well as the data fields for the
     * corresponding samples
     * 
     * @param fileName
     * @param samples
     * @param finalData
     * @throws IOException
     */
    public static void toFile(String fileName, String[] samples, ArrayList<String[]> finalData) throws IOException {
        BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
        out.write("#1.2\n");

        // no. of rows:
        out.write(""+finalData.size());
        out.write('\t');
        out.write(""+samples.length);
        out.write("\nName\tDescription");
        for (String samp: samples){
            out.write('\t');
            out.write(samp);
        }
        out.write('\n'); // header lines done

        for (int row = 0; row < finalData.size(); row++){
            String[] cells = finalData.get(row);
            out.write(cells[0]);
            for (int j = 1; j < cells.length;j++){
                out.write('\t');
                out.write(cells[j]);
            }
            out.write('\n');
        }
        out.close();
    }


    public ArrayList<GCTEntry> getTranscripts() {
        ArrayList<GCTEntry> list = new ArrayList<GCTEntry>();
        for (int i = 0 ; i < nRows; i ++){
            list.add(new GCTEntry(this.ids[i],this.desc[i],this.data[0][i]));
        }
        return list;
    }


    public ArrayList<GCTEntry> getTranscripts(String sampId) {
        ArrayList<GCTEntry> list = new ArrayList<GCTEntry>();

        int index = -1;
        for (int i = 0 ; i < nSamps; i++){
            if (this.samples[i].equals(sampId)){
                index = i;
                break;
            }
        }
        if (index <0){
            throw new RuntimeException("Sample id not found in GCT file: " + sampId+". Samples were: " + Arrays.toString(samples));
        }
        for (int i = 0 ; i < nRows; i ++){
            list.add(new GCTEntry(this.ids[i],this.desc[i],this.data[index][i]));
        }
        return list;
    }


    /**
     * Makes a copy of this GCT file but includes only those rows that have an id matching the modelIds set
     *
     * @param modelIds
     * @return
     */
    public GCTFile filterById(Set<String> modelIds) {
        int newLength = 0;
        for (String id: this.ids){
            if (modelIds.contains(id)){
                newLength++;
            }
        }

        GCTFile newGCT = new GCTFile(newLength , this.nSamps);

        // add ids and desc
        int newRowDex = 0;
        for (int i = 0 ; i < this.nRows; i++){
            if (modelIds.contains(this.ids[i])){
                newGCT.ids[newRowDex] = this.ids[i];
                newGCT.desc[newRowDex] = this.desc[i];
                newRowDex++;
            }
        }
        // add the data
        for (int i = 0;i < this.nSamps; i++){
            // add the header
            newGCT.samples[i]=this.samples[i];
            // add the data by row
            newRowDex=0;
            for (int j = 0; j < nRows; j++){
                if (modelIds.contains(this.ids[j])){
                    newGCT.data[i][newRowDex] = this.data[i][j];
                    newRowDex++;
                }
            }
        }
        return newGCT;
    }


    public int countExpressed(double lowerExprCutoff) {
        int count =0;
        for (int row = 0 ; row < data[0].length; row++){
            if (data[0][row] >= lowerExprCutoff)count++;
        }
        return count;
    }


    public GCTFile combinePreAligned(GCTFile gct) {
        GCTFile newGCT = new GCTFile(this.nRows , this.nSamps+gct.nSamps);

        // add the sample and description fields:
        for (int row = 0 ; row < this.nRows; row++){
            newGCT.ids[row] = this.ids[row];
            newGCT.desc[row] = this.desc[row];
        }

        // add the sample names and data from 'this' objected
        int newSampDex = 0;
        for (;newSampDex < this.nSamps; newSampDex++){
            // add the header
            newGCT.samples[newSampDex]=this.samples[newSampDex];

            // add the data by row
            for (int j = 0; j < nRows; j++){
                newGCT.data[newSampDex][j] = this.data[newSampDex][j];
            }
        }

        // now add the info from the second gct
        int i = 0;
        for (;newSampDex < newGCT.nSamps; newSampDex++){
            //System.out.println("Adding additional col: " + gct.samples[i]);
            // add the header
            newGCT.samples[newSampDex]=gct.samples[i];
            // add the data by row
            for (int j = 0; j < nRows; j++){
                newGCT.data[newSampDex][j] = gct.data[i][j];
            }
            i++;
        }
        return newGCT;
    }


    public void setSampleId(int index, String sampId) {
        this.samples[index] = sampId;
    }


    public void toFile(String outfile, int digits) throws IOException{
        BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
        out.write("#1.2\n");

        // no. of rows:
        out.write(""+this.nRows);
        out.write('\t');
        out.write(""+nSamps);
        out.write("\nName\tDescription");
        for (String samp: samples){
            out.write('\t');
            out.write(samp);
        }
        out.write('\n'); // header lines done

        for (int row = 0; row < nRows; row++){
            out.write(ids[row]);
            out.write('\t');
            out.write(desc[row]);

            for (int j = 0 ; j < nSamps; j++){
                out.write('\t');
                out.write(String.format("%."+digits+"f",data[j][row]));
            }
            out.write('\n');
        }
        out.close();
    }


    public static class GCTEntry {
        String id, desc;
        double val;

        GCTEntry(String id_, String desc_, double v) {
            id = id_;
            desc = desc_;
            val = v;
        }

        public static ArrayList<GCTEntry> sortExpressed(ArrayList<GCTEntry> transcripts, float lowerCutoff) {
            TreeSet<GCTEntry> tree = new TreeSet<GCTEntry>(new GCTEntrySortComparator());
            for (GCTEntry e: transcripts){
                if (e.val >= lowerCutoff){
                    tree.add(e);
                }
            }
            ArrayList<GCTFile.GCTEntry> expr = new ArrayList<GCTFile.GCTEntry>(tree);
            return expr;
        }

        public String getId() {
            return id;
        }
    }


    public static class GCTEntrySortComparator implements Comparator<GCTEntry> {

        public int compare(GCTEntry gctEntry, GCTEntry gctEntry1) {
            if (gctEntry.val > gctEntry1.val) {
                return 1;
            }
            return -1;
        }
    }

    public double[] getColumnData(int sampleIndex) {
        return data[sampleIndex];
    }

    public double[] getColumnData(String sampleName) {
        for (int i = 0 ; i < samples.length; i++) {
            if (samples[i].equals(sampleName)) {
                return data[i];
            }
        }
        return null;
    }

    public void setDesc(HashMap<String,String> idToDescMap) {
        for (int i = 0 ; i < ids.length; i++) {
            String desc = idToDescMap.get(ids[i]);
            if (desc!=null){
                this.desc[i] = desc;
            }
        }
    }

}
