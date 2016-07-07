package org.broadinstitute.cga.tools;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.StringTokenizer;
import java.util.Vector;

public class IOTools {

    /**
     * This simply ensures that the return string ends in a slash,
     * whether or not the supplied string ends in a slash
     *
     * @param _dirpath
     * @return
     */
    static private String parseDirPath(String _dirpath) {
        if (!_dirpath.endsWith("/")) _dirpath+="/";
        return _dirpath;
    }


    /**
     * Deletes a directory recursively with all the files in it
     * 
     * @param dir
     * @return
     */
    public static boolean deleteDirectory(File dir) {
        if (dir.isDirectory()) {
            String[] children = dir.list();
            for (int i=0; i<children.length; i++) {
                boolean success = deleteDirectory(new File(dir, children[i]));
                if (!success) {
                    return false;
                }
            }
        }
        // The directory is now empty so delete it
        return dir.delete();
    } 


    /**
     * Creates a directory
     * @param _DirectoryPath
     * @return
     */
    static public boolean createDirectory(String _DirectoryPath) {
        if (!directoryExists(_DirectoryPath)) {
            File o=new File(parseDirPath(_DirectoryPath));
            return o.mkdir();
        } else {
            return true;
        }
    }


    /**
     * Checks if a directory exists
     * @param _DirectoryPath
     * @return
     */
    static public boolean directoryExists(String _DirectoryPath) {
        _DirectoryPath=parseDirPath(_DirectoryPath);
        File dir=new File(_DirectoryPath);
        return dir.exists();
    }


    /**
     * Gets the directory and returns the files contained on
     * an array
     * @param _DirectoryPath
     * @return
     */
    static public String[] getDirectory (String _DirectoryPath) {
        File dir =new File(_DirectoryPath);
        return dir.list();
    }


    /**
     * Returns the size of a given file to bytes 
     * @param _fileName
     * @return
     */
    static public long getFileSize(String _fileName) {
        File f=new File(_fileName);
        return f.length();
    }
    /**
     * 
     * @param filename
     * @param delims, if null then default Tokenizer delims are used: " \t\n\r\f"
     * @return
     * @throws IOException
     */
    static public ArrayList<String> fileToList(String filename, String delims) throws IOException {

        if (delims == null){
            delims = " \t\n\r\f"; // white space
        }

        ArrayList<String> vec = new ArrayList<String>();
        BufferedReader in = new BufferedReader ( new FileReader(filename));
        String line = in.readLine();

        while (line != null) {
            StringTokenizer toks = new StringTokenizer(line, delims);
            while (toks.hasMoreTokens()) {
                String id =toks.nextToken(); 
                vec.add(id);
            }
            line  = in.readLine();
        }
        in.close();
        return vec;
    }


    static public void urlToFile(String urlStr, String downloadedFileName) throws IOException {
        URL url = new URL(urlStr);
        BufferedInputStream in = new BufferedInputStream(url.openStream());
        BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(downloadedFileName));

        byte[] b = new byte[1000];
        int c = in.read(b);
        while (c >=0){
            out.write(b,0,c);
            c = in.read(b);
        }
        in.close();
        out.close();
    }


	/**
	 * Reads a table from a file into a vector of string arrays.
	 * The delimiter for the fields must be put into the regExp parameter
	 * (utilizes the String.split(regEx) function
	 * 
	 * @param filename
	 * @param regExp
	 * @return
	 * @throws IOException
	 */
    public static ArrayList<String[]> fileToListofArrays(String filename,
            String regExp) throws IOException {

        File f = new File(filename);
        if (!f.exists()){
            return null;
        }

        if (regExp == null){
            regExp = "\\t";
        }

        ArrayList<String[]> vec = new ArrayList<String[]>();
        BufferedReader in = new BufferedReader (new FileReader(f));
        String line = in.readLine();

        while (line != null) {
            String [] split = line.split(regExp,-1);
            vec.add(split);
            line  = in.readLine();
        }
        in.close();
        return vec;
    }


    /**
     * Reads a table from a file into a vector of string arrays.
     * The delimiter for the fields must be put into the regExp parameter
     * (utilizes the String.split(regEx) function
     *
     * @param inStream
     * @param regExp
     * @return
     * @throws IOException
     */
    public static ArrayList<String[]> fileToListofArrays(InputStream inStream,
                                                         String regExp) throws IOException {

        if (regExp == null){
            regExp = "\\t";
        }

        ArrayList<String[]> vec = new ArrayList<String[]>();
        BufferedReader in = new BufferedReader ( new InputStreamReader(inStream));
        String line = in.readLine();

        while (line != null){
            String [] split = line.split(regExp,-1);
            vec.add(split);
            line  = in.readLine();
        }
        in.close();
        return vec;
    }


    public static String getString(InputStream stream) throws IOException {
        return getString(new InputStreamReader(stream));
    }


    public static String getString(Reader reader) throws IOException {
        BufferedReader in = new BufferedReader ( reader);

        String line = in.readLine();
        StringBuilder str = new StringBuilder();
        while (line != null){
            str.append(line);
            str.append("\n");
            line  = in.readLine();
        }
        in.close();
        return str.toString();
    }


    static public String fileToString(String filename) throws IOException {
        return getString(new FileReader (filename));
    }


    static public String fileToString(File file) throws IOException {
        return getString(new FileReader (file));
    }


    public static void stringToFile(String filename, String value) throws IOException {
        BufferedWriter out = new  BufferedWriter(new FileWriter(filename));
        out.write(value);
        out.close();
    }


    public static void listArrayToFile(ArrayList<String[]> rows, String filename) throws IOException{
        BufferedWriter out = new BufferedWriter( new FileWriter(filename));
        for (String[] line: rows) {
            for (String cell: line) {
                out.write(cell);
                out.write('\t');
            }
            out.write('\n');
        }
        out.close();
    }


    /**
     *  Copies src file to dst file.  If the dst file does not exist, it is created
     *
     * @param src
     * @param dst
     * @throws IOException
     */
    public static void copy(File src, File dst) throws IOException {
        if (src.getPath().equals(dst.getPath())) return ; // don't copy if they are the same file

        InputStream in = new FileInputStream(src);
        OutputStream out = new FileOutputStream(dst);

        // Transfer bytes from in to out
        byte[] buf = new byte[1024];
        int len;
        while ((len = in.read(buf)) > 0) {
            out.write(buf, 0, len);
        }
        in.close();
        out.close();
    }


    public static String tableFileToHTML(String docTableFile) throws IOException {
        ArrayList<String[]> table = fileToListofArrays(docTableFile, null);
        if (table == null || table.size() == 0) return null;

        StringBuilder str = new StringBuilder();
        str.append("<table>");
        for (String header: table.get(0)) {
            str.append("<th>").append(header).append("</th>");
        }
        str.append('\n');

        for (int i = 1; i < table.size(); i++) {
            String[] cells = table.get(i);
            str.append("<tr>");
            for (String c: cells) {
                boolean isNumeric = true;
                try {
                    Double.valueOf(c);
                } catch(NumberFormatException e) {
                    isNumeric = false;
                }
                if (isNumeric) {
                    str.append("<td align='right'>").append(c).append("</td>");
                } else {
                    str.append("<td>").append(c).append("</td>");
                }
            }
            str.append("</tr>\n");
        }
        str.append("</table>\n");
        return  str.toString();
    }


    public static void collectionToFile(String filename, Iterable c) throws IOException {
        BufferedWriter out = new BufferedWriter(new FileWriter(filename));
        for (Object item: c) {
            out.write(item.toString());
            out.write('\n');
        }
        out.close();
    }


    /**
     *
     * @param dirName
     * @return true if a directory was created
     */
    public static boolean mkDir(String dirName) {
        File dir = new File(dirName);
        return dir.mkdir();
    }


    public static void gobble(ProcessBuilder pb) throws InterruptedException, IOException {
        Process p = pb.start();
        IOTools.StdOGobbler gobbler = new IOTools.StdOGobbler(p);
        IOTools.StdErrGobbler errGobbler = new IOTools.StdErrGobbler(p);
        gobbler.start();
        errGobbler.start();

        gobbler.join();
        errGobbler.join();
        p.waitFor();
    }


    public static class StdOGobbler extends Thread {
        Process p;
        String filename = null;

        public StdOGobbler (Process pro) {
            super();
            p = pro;
        }

        public void setToFile(String filename){
            this.filename = filename;
        }
        // This method is called when the thread runs
        public void run() {
            try {
                InputStream in = this.selectStream();
                if (filename!= null) {
                    BufferedOutputStream out = new BufferedOutputStream (new FileOutputStream(filename));
                    gobble(in, out);
                    out.close();
                } else {
                    gobble(in, System.out);
                }
                in.close();
            } catch (IOException e) {
                e.printStackTrace();
                throw new RuntimeException ("IO problem gobbling stdout from process (GNUPlot) ", e);
            }
        }

       protected InputStream selectStream() {
           return p.getInputStream();
       }

       private void gobble(InputStream in, OutputStream out) throws IOException {
            int c;
            byte[] bytes = new byte[1000];
            while ((c = in.read(bytes))>0) {
                out.write(bytes,0,c);
            }
        }
    }


    public static class StdErrGobbler extends StdOGobbler {
        public StdErrGobbler (Process pro) {
            super(pro);
        }

        protected InputStream selectStream() {
            return super.p.getErrorStream();
        }
    }

}
