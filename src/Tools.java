package org.broadinstitute.cga.tools;


import javax.mail.Session;
import javax.mail.MessagingException;
import javax.mail.internet.MimeMessage;
import javax.mail.internet.InternetAddress;

import java.util.*;
import java.net.URLEncoder;
import java.net.URLDecoder;
import java.io.*;


/**
 * Created for Hanover Medical School.
 * User: dave
 * Date: 13.05.2005
 * Time: 15:57:13

 This is a set of static methods which perform random, helpful tasks ...


    <pre>
        <Tools>
            <Properties>
                <TEMP type="relativePath">./resources/temp</TEMP>
            </Properties>
        </Tools>
    </pre>
 */

public class Tools {

	/**
	 * This translate a comma separated list of numbers into a vector of numbers
	 * For example:
	 *
	 * 1,2,3-6
	 *
	 * would then contain a vector having all 6 positions
	 *
	 * @param numList  a list of number separated by dashes and or commas
	 * @return
	 */
	static public Vector<Integer> parseNumList(String numList) {
		Vector<Integer> vec = new Vector<Integer>(5);
		StringTokenizer tok = new StringTokenizer(numList, ", ");
		String num;
		while (tok.hasMoreTokens()){
			num = tok.nextToken();
			if (num.lastIndexOf("-") != -1){
				String[] strs = num.split("-");
				if (strs.length != 2) throw new RuntimeException("Could not parse the list of numbers");
				int start = Integer.parseInt(strs[0]);
				int stop = Integer.parseInt(strs[1]);
				for (int i = start; i <= stop; i++) vec.add(new Integer(i));
			}else vec.add(new Integer(num.trim()));
		}
		return vec;
	}

	/**
	 * Takes an array of positions and condenses them using '-' when possible
	 *
	 * @param positions
	 * @return
	 */
	static public String condensePositions (String[] positions) {

		if (positions.length == 0)
			return "";

		int next, cCount = 0;
		boolean consecutive;

		String posList = positions[0];
		for (int i=1; i<positions.length; i++) {
			next = i + (i == positions.length - 1 ? 0 : 1);
			consecutive = Integer.parseInt(positions[i]) - Integer.parseInt(positions[i-1]) == 1;
			if (consecutive) {
				cCount += 1;
				if (Integer.parseInt(positions[next]) - Integer.parseInt(positions[i]) != 1)
					posList = posList + (cCount > 1 ? "-" : ",") + positions[i];
			} else {
				cCount = 0;
				posList = posList + "," + positions[i];
			}
		}

		return posList;
	}

	/**                                                                             
	 * Returns the working directory as given in the system property
	 *
	 * does:  return System.getProperty("user.dir").replace('\\', '/');
	 *
	 * @return
	 */
	static public String getWorkingDirectory() {
		return System.getProperty("user.dir").replace('\\', '/');        
	}

	/**
	 * Converts a collection of objects to string by calling the toString()
	 * method of each object.
	 *
	 * The result has the form: [ hi, bye ]
	 *
	 * @param collection   this is turned into a string
	 * @return the string
	 */
	public static String collectionToString(Collection<? extends Object> collection) {
		if(collection == null || collection.size() == 0) return "[ ]";
		StringBuilder str = new StringBuilder(20);
		str.append('[');
		for (Object o: collection) {
			str.append(o.toString());
			str.append(", ");
		}
		str.delete(str.length()-2,str.length());
		str.append(']');
		return str.toString();
	}

	/**
	 * Converts a collection of objects to string by calling the toString()
	 * method of each object.
	 *	
	 * @param collection   this is turned into a string
	 * @return the string
	 */
	public static String collectionToString(Collection<? extends Object> collection, String delim) {
		if(collection == null || collection.size() == 0) return "";
		StringBuilder str = new StringBuilder(20);
		int last = collection.size();
		int count = 0;
		for (Object o: collection) {
			count++;
			str.append(o.toString());
			if (count < last) str.append(delim);
		}
	
		return str.toString();
	}


	static String ENCODING ="UTF-8";

	/**
	 * Translates a string into application/x-www-form-urlencoded format using a
	 * UTF-8.
	 *
	 * This method uses the supplied encoding scheme to obtain the bytes for unsafe
	 * characters
	 *
	 * @param str
	 * @return
	 * @throws UnsupportedEncodingException
	 */
	static public String encode(String str) throws UnsupportedEncodingException {
		return URLEncoder.encode(str, ENCODING);
	}
	
	/**
	 * Decodes a application/x-www-form-urlencoded string using UFT-8.
	 *
	 *
	 * @param str
	 * @return
	 * @throws UnsupportedEncodingException
	 */
	static public String decode(String str) throws UnsupportedEncodingException {
		if (str == null) return null;
		return URLDecoder.decode(str, ENCODING);
	}

	
	/**
	 * one way encoding to remove special characters (replaced by underscores)
	 * @param str
	 * @return
	 * @throws UnsupportedEncodingException
	 */
	static public String encodeHomeBrew(String str) throws UnsupportedEncodingException {
		str = str.replace('|', '_');
		str = str.replace('\t', '_');
		str = str.replace(' ', '_');
		//TODO add additional replacements
		return str;
	}
	


	
	/**
	 * No idea what this does, or where it's used
	 *
	 * Nektarios????
	 *
	 * @param _toParse
	 * @param _findString
	 * @param delimeter
	 * @return
	 */
	static public String parseString(String _toParse,String _findString,String delimeter)
	{

		String[] split=_toParse.split(_findString);
		int splitPos=0;
		if (split.length==2)
		{
			for (int i=0;i<split[1].length();i++)
			{
				splitPos+=1;					
				if (split[1].charAt(i)==delimeter.charAt(0))
				{
					break;
				}				
			}
			return split[1].substring(0,splitPos-1);
		}
		else
		{
			return null;	
		}				
	}

	/**
	 * Performs the initializations necessary when an application starts.
	 *
	 * 1) Initializes the StaticProperties class using the conf/config_order.txt file
	 * 2) Initializes the ATMLogger using log4j xml file, as defined in static properties
	 *
	 */
	static public void initialize(){
		try{
			if (StaticProperties.properties == null) {
				StaticProperties.init();

				ATMLogger.init(StaticProperties.get("LOG4J_CONF_FILE", ATMLogger.class),
						StaticProperties.ROOT_DIR );
				ATMLogger.debug(StaticProperties.showSettings());
				ATMLogger.debug("Working directory: " + Tools.getWorkingDirectory());
				ATMLogger.debug("ATM_ROOT: " + StaticProperties.ROOT_DIR);
			}
		}catch (Exception e){
			e.printStackTrace();
			throw new RuntimeException (e);
		}
	}


	/**
	 * initializes system and forces reloading of application settings as given
	 *
	 * @param forceReload      forces the reloading of the static properties if true
	 */
	static public void initialize(boolean forceReload){
		if (forceReload) {
			StaticProperties.properties = null;
		}
		initialize();
	}

	/**
	 * Removes any white space anywhere within this string
	 *
	 * White space is determined as having a greater ascii value than ' ' (blank)
	 *
	 * @param s            a string
	 * @return
	 */
	static public String removeWhiteSpace(String s) {
		int l = s.length();
		char c;
		StringBuilder str = new StringBuilder(l);
		for (int i = 0 ; i <l; i++){
			c = s.charAt(i);
			if (c > ' ')
				str.append(c);
		}
		return str.toString();
	}

	/**
	 * returns the path to the TEMP directory as defined in config.xml
	 *
	 * @return
	 */
	static public String getTempDir() {
		return StaticProperties.get("TEMP", Tools.class);
	}

	/**
	 * Send an email. Utilizes the SMTP server given in config.xml
	 * 
	 * @param senderEmail          from
	 * @param recipient            to
	 * @param subject              subject
	 * @param body                 body
	 * @throws MessagingException
	 */
	static public void sendEmail(String senderEmail, String recipient, String subject,
			String body) throws MessagingException{

		Properties props = new Properties();
		props.put("mail.smtp.host", StaticProperties.get("SMTP_HOST",Tools.class));
		Session session = Session.getDefaultInstance(props, null);
		MimeMessage message = new MimeMessage(session);
		message.setFrom( new InternetAddress(senderEmail));
		message.addRecipient(MimeMessage.RecipientType.TO, new InternetAddress(recipient));
		message.setSubject(subject);
		message.setText(body);
		javax.mail.Transport.send(message);

	}


	/**
	 * removes fasta header, if present
	 * removes all white space
	 * returns sequence in uppercase letters
	 * @param text
	 * @return
	 */
	static public String textAreaToSequence(String text){
		text = text.trim();
		if (text.startsWith(">")){
			text = text.substring(text.indexOf('\n')+1);
		}
		text = text.replace("\t", "").replace(" ", "").replace("\r", "").replace("\n", "");
		return text.toUpperCase();
	}

	/**
	 * 
	 * @param header
	 * @param sequence
	 * @param lineLength  zero o less to use the default line length (60)
	 * @return
	 */
	static public String toFasta(String header, String sequence, int lineLength){
		if (lineLength <1){
			lineLength = 60;
		}

		StringBuilder str = new StringBuilder();
		if (header !=null){
			str.append('>').append(header ).append('\n');
		}
		int i = 0;
		int l = sequence.length();
		do {
			int end = i+lineLength;
			if ( i+lineLength  > l){
				end = l;
			}
			str.append(sequence.substring(i,end)).append('\n');
			i = end;
		} while (i < l);
		str.append( "\n");
		return str.toString();

	}


	public static String getStackTrace(Throwable aThrowable) {
		final Writer result = new StringWriter();
		final PrintWriter printWriter = new PrintWriter(result);
		aThrowable.printStackTrace(printWriter);
		return result.toString();
	}

    public static String stripFileName(String filePath) {
        File f = new File(filePath);
        return f.getName();
    }


    /**
     * Split the string into tokens separated by the given delimiter.  Profiling has
     * revealed that the standard string.split() method typically takes > 1/2
     * the total time when used for parsing ascii files.
     * @author based on code written by Jim Robinson
     * @param aString
     * @param delim
     * @return
     */
    public static ArrayList splitString(String aString, String delim)
    {
        aString = aString.trim();
        int start = 0;
        int end = aString.indexOf(delim);
        ArrayList<String> tokens = new ArrayList();

        while (end > 0) {
            tokens.add(aString.substring(start, end));
            start = end + 1;
            end = aString.indexOf(delim, start);
        }

        tokens.add(aString.substring(start));

        return tokens;
    }
}


