package org.broadinstitute.cga.tools;

import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.zip.GZIPInputStream;

/**
 * User: DeLuca
 * Date: 21.02.2006
 * Time: 13:05:04
 * This code will blow your mind!
 */
public class DynamicStreams {
    /**
     * When s begings with "http://" or "ftp://" then it returns an InputStream from a URL
     *     otherwise the source is a File
     * Additionally:
     *     s ends with ".gz" then open from GZIP Stream
     *	s ends with ".zip" then open from ZIP stream
     *	WHAT HAPPENS WHEN THERE ARE DIRECTORIES AND MULTIPLE FILES IN ZIP ARCHIVE?????
     * @param s     filename, or url
     * @return      The resulting stream is buffered
     * @throws java.net.MalformedURLException
     * @throws java.io.IOException
     */	
    static public InputStream getDynamicInputStream(String s) throws MalformedURLException, IOException {
        BufferedInputStream in = null;
        String lower  = s.trim().toLowerCase();

        if (lower.startsWith("http://")) {
            in = new BufferedInputStream((new URL(s)).openStream());
        } else {
            in = new BufferedInputStream(new FileInputStream(s));
        }

        if (lower.endsWith(".gz")) {
            return new GZIPInputStream(in);
        }
        return in; // iw was a non-compressed file
    }

    /**
     * Wrapps a reader around the stream returned by DynamicStreams.getDynamicInputStream
     * @param s
     * @return
     * @throws IOException
     */
    static public BufferedReader getDynamicReader(String s) throws IOException {
        return new BufferedReader(new InputStreamReader(getDynamicInputStream(s)));
    }

    static public BufferedWriter getBufferedWriter(String s) throws IOException {
        return new BufferedWriter(new FileWriter(s));
    }

    /**
     * reads the rest of the stream so that the process can end
     * @param in
     */
    static public  void eatStream(Reader in) throws IOException {
        while (in.read() != -1);
    }

    /**
     * eats the number of ocurrences of '\n' specified by n
     * @param n
     * @param in
     * @return returns true if the end of file has been reached
     * @throws IOException
     */
    static public boolean eatLines(int n, Reader in) throws IOException {
        int lineCount = 0;
        int c;
        do {
            c = in.read();
            if (((char)c) == '\n') {
                lineCount ++;
            }
        } while (c != -1 && lineCount < n);
        return c == -1;
    }

    /**
     * reads from the Reader until the character c has been found
     * @param c
     * @param in
     * @return  returns true if the end of the file has been reached
     * @throws IOException
     */
    static public boolean advanceTo(char c, Reader in) throws IOException {
        int intC = (int) c;
        int thisC;
        do {
            thisC = in.read();
        } while (thisC != -1 && thisC != intC);
        return thisC == -1;
    }

    public static final char EOF = (char) 4;

    public static void main(String [] args) {
        try {
            Tools.initialize();
            System.out.println("Getting stream");
            InputStream in = DynamicStreams.getDynamicInputStream("ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/Alignments_Rel_2.16.ZIP");
            System.out.println("Writing to file");
            FileOutputStream out = new FileOutputStream("f:/align.zip");
            int c = in.read();
            if (c == -1) {
                System.out.println("Error, the file is empty before we started!!!!");
            }
            while (c != -1) {
                out.write(c);
                c = in.read();
            }
            System.out.println("\nClosing");
            out.close();
            in.close();
            System.out.println("Done.");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * does not close either stream when finished
     * @param from
     * @param to
     * @throws IOException
     */
    public static void redirect(InputStream from, OutputStream to) throws IOException {
        byte[] buff = new byte[1000];
        int n = from.read(buff);
        while (n >=0) {
            to.write(buff,0,n);
            n = from.read(buff);
        }
    }

}
