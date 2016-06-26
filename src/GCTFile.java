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

    public GCTFile(String filename) throws IOException{
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
        for (int i = 2; i < headers.length;i++){
            samples[i-2] = headers[i];
        }

        line = in.readLine(); // actual first data line


        int row = 0;
        while (line!= null){
            String [] split = line.split("\\t"); // form: Accname  desc    sampleName
            ids[row] = split[0];
            desc[row] = split[1];

//            System.out.println("id\t"+ids[row]);
//            System.out.println("desc\t"+desc[row]);
//            System.out.println("split\t" + Arrays.toString(split));
            for (int i = 2; i < split.length; i++){
                float val = Float.valueOf(split[i]);
//                System.out.println("val:\t"+val);
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
        while (table.get(start)[0].startsWith("#")){
            start++;
        }
        String[] header = table.get(start);
        start++;

        nRows = table.size() - start;
        nSamps = header.length -1 ;

        ids = new String[nRows];
        desc = new String[nRows];
        samples = new String[nSamps];
        data = new double[nSamps][nRows];

        for (int i = 1; i < header.length;i++){
            samples[i-1] = header[i];
        }


        int row = 0;
        for (int oldRow = start; oldRow < table.size(); oldRow++){
            String [] split = table.get(oldRow);
            ids[row] = split[0];
            desc[row] = "";

            for (int i = 1; i < split.length; i++){
                float val = Float.valueOf(split[i]);
//                System.out.println("val:\t"+val);
                data[i-1][row]= val;

            }
            row++;
        }


    }

    /**
     * This parsese the GCT file. The colapsing happens just throught he nature of the returned structure being
     * a map, which can only contain one value per key
     *
     * @param filename
     * @return acc to data map
     * @throws IOException
     */
    private static HashMap<String, Double> getCollapsedData(String filename)throws IOException {
        BufferedReader in = new BufferedReader (new FileReader(filename));
        HashMap<String,Double> data = new HashMap<String, Double>();

        Hashtable<String,Double> table = new Hashtable<String,Double>();
        String line = in.readLine(); // that weird #1.2 line
        line = in.readLine(); // line with number of rows and columns
        line = in.readLine(); // headers
        line = in.readLine(); // actual first data line


        while (line!= null){
            String [] split = line.split("\\t"); // form: Accname  desc    sampleName
            String acc = split[0];
            double value = Double.parseDouble(split[2]);

            data.put(acc,value);

            line = in.readLine();
        }


        in.close();

        return data;
    }


    private static void alignExpressionData(HashMap<String, Double> expValues, String refGCT) throws IOException{


        BufferedWriter out = new BufferedWriter(new FileWriter("/Users/ddeluca/Documents/GTEx/Factorial/GCT/K562_aligned.txt"));

        BufferedReader in = new BufferedReader (new FileReader(refGCT));

        in.readLine();
        in.readLine();
        in.readLine();
        String line = in.readLine();

        while(line!=null){

            String gene = line.split("\\t")[1];
            if (expValues.containsKey(gene)){
                out.write(expValues.get(gene).toString()); out.write('\n');
            }else {
                out.write(gene); out.write("\tNA\n");
            }

            line = in.readLine();
        }
        out.close();
        in.close();


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
    public GCTFile combine(GCTFile gct, boolean useDescFromPassedObj, boolean useDescFromThisObj)  {
        GCTFile newGCT = new GCTFile(this.nRows , this.nSamps+gct.nSamps);

        int newSampDex = 0;

        for (;newSampDex < this.nSamps; newSampDex++){
            // add the header
            newGCT.samples[newSampDex]=this.samples[newSampDex];

            // add the data by row
            for (int j = 0; j < nRows; j++){
                newGCT.data[newSampDex][j] = this.data[newSampDex][j];
                newGCT.ids[j] = this.ids[j];
                newGCT.desc[j] = this.desc[j];
            }
        }

        // now add the info from the second gct
        HashMap<String,double[]> indexed = null;
        if(useDescFromPassedObj)
        {
            indexed = gct.getIndexedByDesc();
            int i = 0;
            for (;newSampDex < newGCT.nSamps; newSampDex++){
                //System.out.println("Adding additional col: " + gct.samples[i]);
                // add the header
                newGCT.samples[newSampDex]=gct.samples[i];
                // add the data by row
                for (int j = 0; j < nRows; j++){
                    double[] row;
                    if(useDescFromThisObj)
                    {
                        row = indexed.get(newGCT.desc[j]);
                    }
                    else
                    {
                        row = indexed.get(newGCT.ids[j]);
                    }
                    if (row ==null){
                        newGCT.data[newSampDex][j] = Float.NaN;
                    }else{
                          newGCT.data[newSampDex][j] = row[i];
                    }
                }
                i++;
            }
        }
        else
        {
            indexed = gct.getIndexedById();
            int i = 0;
            for (;newSampDex < newGCT.nSamps; newSampDex++){
                //System.out.println("Adding additional col: " + gct.samples[i]);
                // add the header
                newGCT.samples[newSampDex]=gct.samples[i];
                // add the data by row
                for (int j = 0; j < nRows; j++){
                    double[] row;
                    if(useDescFromThisObj)
                    {
                        row = indexed.get(newGCT.desc[j]);
                    }
                    else
                    {
                        row = indexed.get(newGCT.ids[j]);
                    }
                    if (row ==null){
                        newGCT.data[newSampDex][j] = Float.NaN;
                    }else{
                         newGCT.data[newSampDex][j] = row[i];
                    }
                }
                i++;
            }
        }


        return newGCT;


    }


    /**
     * gct is the rnaseq file with ensemblids
     * @param gct
     * @return
     * @throws IOException
     */
    public GCTFile combineSpecial(GCTFile gct) throws IOException {
        // first we have to trim off the version numbers:

        for(int i = 0 ; i < gct.ids.length; i++){
            gct.ids[i] = gct.ids[i].substring(0,gct.ids[i].indexOf('.'));
            //System.out.println(gct.ids[i]);
        }

        // create a new GCT file to contain the combined data

        GCTFile newGCT = new GCTFile(this.nRows , this.nSamps+gct.nSamps);

        int newSampDex = 0;

        for (;newSampDex < this.nSamps; newSampDex++){
            // add the header
            newGCT.samples[newSampDex]=this.samples[newSampDex];

            // add the  data from 'this' object, row by row
            for (int j = 0; j < nRows; j++){
                newGCT.data[newSampDex][j] = this.data[newSampDex][j];
                newGCT.ids[j] = this.ids[j]; // ok, this is silly to do over and over ...
                newGCT.desc[j] = this.desc[j];
            }
        }

        // now add the info from the second gct
        HashMap<String,double[]> indexed = null;
        indexed = gct.getIndexedById();
        System.out.println("Indexed rows: " + indexed.size());
        for (String id: indexed.keySet()){
            // id is ensemble id from the rna-seq. now we need to find the row that contains this id in its description

            int matchedRow = -1;

            for (int j = 0 ; j < nRows; j++){
                if (desc[j].contains(id)){
                    matchedRow = j;
                    break;
                }
            }
            if (matchedRow>=0){
                int i =0;
                for (int col = newSampDex
                             ;col < newGCT.nSamps; col++){
                    // add the data by row
                    double[] row = indexed.get(id);
                    if (row ==null){
                        newGCT.data[col][matchedRow] = Float.NaN;
                    }else{
                        newGCT.data[col][matchedRow] = row[i];
                    }
                    i++;
                }
            }
        }

        int i = 0;
        for (int col = newSampDex;col < newGCT.nSamps; col++){
            // add the header
            newGCT.samples[col]=gct.samples[i];
            i++;

        }


        return newGCT;


    }


    /**
     * Returns a new GCT which has the sample and values of this object, but the ordering and gene list
     * of the passed paramater object. The matching occurs via the description field (gene).
     *
     *  If there is no equivalent row in this GCTFile, then a row of Float.NaN is added
     *
     * @param gct
     * @param keepMissingRows In the case that there are rows in gct but not in 'this', then this determines whether
     *              to keep add empty NaN rows or whether to skip the row
     * @return
     * @throws IOException
     */
    public GCTFile sortToMatchById(GCTFile gct, boolean keepMissingRows)  {

        HashMap<String,double[]> indexed = this.getIndexedById();
        System.out.println("Contains UNC5A? " + indexed.containsKey("UNC5A"));

//        try{
//        IOTools.collectionToFile("/Users/ddeluca/Documents/GTEx/AffyCombined/debug"+new Date().getTime()+".txt",indexed.keySet());
//        }catch(Exception e ) {e.printStackTrace();}
        int newRowCount = gct.nRows;
        if (!keepMissingRows){
            newRowCount = 0;
            for (int j = 0; j < gct.nRows; j++){
                if (indexed.containsKey(gct.ids[j])){
                    newRowCount++;
                }else{
                    //System.out.println("not found:\t" +gct.ids[j]);
                }
            }
        }

        System.out.println("Row count reduced from " + gct.nRows +" to " + newRowCount);

        GCTFile newGCT = new GCTFile(newRowCount , this.nSamps);

        // add samples

        for (int newSampDex = 0;newSampDex < this.nSamps; newSampDex++){
            // add the header
            newGCT.samples[newSampDex]=this.samples[newSampDex];
        }

        // add gene and descriptors from passed gct
        int newRowDex=0;
        for (int j = 0; j < gct.nRows; j++){
            if (keepMissingRows || indexed.get(gct.ids[j])!=null){
                newGCT.ids[newRowDex] = gct.ids[j];
                newGCT.desc[newRowDex] = gct.desc[j];
                newRowDex++;
            }
        }


        int i = 0;

        for (int newSampDex = 0;newSampDex < newGCT.nSamps; newSampDex++){
            //System.out.println("Adding additional col: " + gct.samples[i]);

            for (int j = 0; j < newGCT.nRows; j++){

                double[] row = indexed.get(newGCT.ids[j]);
                if (row !=null){
                    newGCT.data[newSampDex][j] = indexed.get(newGCT.ids[j])[i];
                }else {
                    newGCT.data[newSampDex][j] = Float.NaN;
                    if(!keepMissingRows){
                        System.out.println("Error: missing rows found when keepMissingRows set to false");
                    }
                }


            }
            i++;
        }

        return newGCT;


    }



    /**
     * Returns a new GCT which has the sample and values of this object, but the ordering and gene list
     * of the passed paramater object. The matching occurs via the description field (gene).
     *
     *  If there is no equivalent row in this GCTFile, then a row of Float.Nan
     * @param gct
     * @return
     * @throws IOException
     */
    public GCTFile sortToMatch(GCTFile gct, boolean useDesc)  {
        GCTFile newGCT = new GCTFile(gct.nRows , this.nSamps);

        // add samples

        for (int newSampDex = 0;newSampDex < this.nSamps; newSampDex++){
            // add the header
            newGCT.samples[newSampDex]=this.samples[newSampDex];
        }

        // add gene and descriptors from passed gct
        for (int j = 0; j < gct.nRows; j++){
            newGCT.ids[j] = gct.ids[j];
            newGCT.desc[j] = gct.desc[j];
        }


        // now add data fromt this gct in the order of the passed gct
        if(useDesc)
        {
            HashMap<String,double[]> indexed = this.getIndexedByDesc();

            int i = 0;

            for (int newSampDex = 0;newSampDex < newGCT.nSamps; newSampDex++){
                //System.out.println("Adding additional col: " + gct.samples[i]);

                // add the data by row
                for (int j = 0; j < gct.nRows; j++){
                    double[] row = indexed.get(gct.desc[j]);
                    if (row ==null){
                        newGCT.data[newSampDex][j] = Float.NaN;
                    }else{
                        newGCT.data[newSampDex][j] = indexed.get(newGCT.desc[j])[i];
                    }


                }
                i++;
            }
        }
        else
        {
            HashMap<String,double[]> indexed = this.getIndexedById();

            int i = 0;

            for (int newSampDex = 0;newSampDex < newGCT.nSamps; newSampDex++){
                //System.out.println("Adding additional col: " + gct.samples[i]);

                // add the data by row
                for (int j = 0; j < gct.nRows; j++){
                    double[] row = indexed.get(gct.ids[j]);
                    if (row ==null){
                        newGCT.data[newSampDex][j] = Float.NaN;
                    }else{
                        newGCT.data[newSampDex][j] = indexed.get(newGCT.ids[j])[i];
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

            for (int j = 0 ; j < nSamps; j++){
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
        for (int i = 0 ; i< nRows; i++){
            double[] row = new double[nSamps];

            for (int j = 0 ; j < nSamps; j++){
                row[j]=data[j][i];
            }
            if (ids[i] == null){
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
            //System.out.println("ids length: " + ids.length + " first val: " +ids[row]);
            out.write(ids[row]);
            out.write('\t');
            out.write(desc[row]);

            for (int j = 0 ; j < nSamps; j++){

//                System.out.println("sampdex:\t"+j);
//                System.out.println("rowdex:\t"+row);
//                System.out.println("rows for samp: " + data[j].length);
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

    public GCTFile collapse() {
        HashMap<String,ArrayList<MathSet>> descriptions = new HashMap<String,ArrayList<MathSet>>();
        for (int rowDex = 0 ; rowDex < this.nRows; rowDex++){
            if (desc[rowDex] != null && !desc[rowDex].replace("-","").equals("")){
                ArrayList<MathSet> mathCols = descriptions.get(desc[rowDex]);
                if (mathCols == null){
                    mathCols = new ArrayList<MathSet>();
                    for (int colDex = 0 ; colDex < nSamps; colDex++){
                        mathCols.add(new MathSet());
                    }
                    descriptions.put(desc[rowDex],mathCols);
                }
                for (int colDex = 0 ; colDex < nSamps; colDex++){
                    mathCols.get(colDex).add((float)this.data[colDex][rowDex]);
                }
            }
        }

        //make new GCTFile based on the descriptions object
        GCTFile newGCT = new GCTFile(descriptions.size(),nSamps);

        /* VALUES TO SET:  String[] ids ;  String desc[] ;String[] samples;  double [][] data; // [samps,rows]
        */

        for (int i = 0; i < newGCT.nSamps; i++){
            newGCT.samples[i] = this.samples[i];
        }

        int row = 0;
        for (String d: descriptions.keySet()){
            if (d == null) {
                System.out.println("How is this null?????? after collapse");
            }
            newGCT.ids[row] = d;
            newGCT.desc[row] = d;

            ArrayList<MathSet> mathRow = descriptions.get(d);
            for (int col = 0; col < nSamps; col++){
                newGCT.data[col][row] = mathRow.get(col).getMean();
            }
            row++;
        }

        if (row == 0){
            System.out.println("Warning: collapsed to zero rows.");
        }
        return newGCT;

    }

    public GCTFile collapseByMax() {
        HashMap<String,ArrayList<MathSet>> descriptions = new HashMap<String,ArrayList<MathSet>>();
        for (int rowDex = 0 ; rowDex < this.nRows; rowDex++){
            if (desc[rowDex] != null && !desc[rowDex].replace("-","").equals("")){
                ArrayList<MathSet> mathCols = descriptions.get(desc[rowDex]);
                if (mathCols == null){
                    mathCols = new ArrayList<MathSet>();
                    for (int colDex = 0 ; colDex < nSamps; colDex++){
                        mathCols.add(new MathSet());
                    }
                    descriptions.put(desc[rowDex],mathCols);
                }
                for (int colDex = 0 ; colDex < nSamps; colDex++){
                    mathCols.get(colDex).add((float)this.data[colDex][rowDex]);
                }
            }
        }

        //make new GCTFile based on the descriptions object
        GCTFile newGCT = new GCTFile(descriptions.size(),nSamps);

        /* VALUES TO SET:  String[] ids ;  String desc[] ;String[] samples;  double [][] data; // [samps,rows]
        */

        for (int i = 0; i < newGCT.nSamps; i++){
            newGCT.samples[i] = this.samples[i];
        }

        int row = 0;
        for (String d: descriptions.keySet()){
            if (d == null) {
                System.out.println("How is this null?????? after collapse");
            }
            newGCT.ids[row] = d;
            newGCT.desc[row] = d;

            ArrayList<MathSet> mathRow = descriptions.get(d);
            for (int col = 0; col < nSamps; col++){
                newGCT.data[col][row] = mathRow.get(col).getMax();
            }
            row++;
        }

        if (row == 0){
            System.out.println("Warning: collapsed to zero rows.");
        }
        return newGCT;

    }


    public int getNumberOfRows() {
        return nRows;
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
            //System.out.println("ids length: " + ids.length + " first val: " +ids[row]);
            out.write(ids[row]);
            out.write('\t');
            out.write(desc[row]);

            for (int j = 0 ; j < nSamps; j++){

//                System.out.println("sampdex:\t"+j);
//                System.out.println("rowdex:\t"+row);
//                System.out.println("rows for samp: " + data[j].length);
                out.write('\t');
                out.write(String.format("%."+digits+"f",data[j][row]));
            }
            out.write('\n');
        }
        out.close();

    }

    /**
     * removes sample columnes to contain only ones that are present in the samples array
     * @param samples if a sample is present in this array AND in this object then it will be retained
     */
    public void subSamples(String[] samples) {

        ArrayList<Integer> keepIndex = new ArrayList<Integer>();
        for (int i = 0; i <this.samples.length;i++){
            // check if this samp is present in the array:
            for (int j = 0; j< samples.length; j++){
                if (this.samples[i].equals(samples[j])){
                    keepIndex.add(i);
                    break;
                }
            }
        }

        // now we know what to keep:

        String[] newSamps = new String[keepIndex.size()];
        double[][] newData = new double[keepIndex.size()][];

        int i = 0;
        for (int iKeep: keepIndex){
            newSamps[i] = this.samples[iKeep];
            newData[i] = this.data[iKeep];
            i++;
        }

        this.nSamps = newSamps.length;
        this.samples = newSamps;
        this.data = newData;
    }

    public String[] getIds() {
        return ids;
    }


    public static class GCTEntry{
        String id, desc;
        double val;

        GCTEntry(String id_, String desc_, double v){
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
            if (gctEntry.val > gctEntry1.val){
                return 1;
            }
            return -1;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }

    public double[] getColumnData(int sampleIndex){
        return data[sampleIndex];
    }

    public double[] getColumnData(String sampleName){
        for (int i = 0 ; i < samples.length; i++){
            if (samples[i].equals(sampleName)){
                return data[i];
            }
        }
        return null;
    }

    public void setDesc(HashMap<String,String> idToDescMap){
        for (int i = 0 ; i < ids.length; i++){
            String desc = idToDescMap.get(ids[i]);
            if (desc!=null){
                this.desc[i] = desc;
            }
        }
    }

}
