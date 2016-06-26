package org.broadinstitute.cga.tools;

import org.jdom.Attribute;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;

import java.util.*;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.BufferedReader;

/**
 * User: DeLuca
 * Date: 17.01.2006
 * Time: 08:27:17
 *

 *
 */
public class StaticProperties {
	protected static String ROOT_ENV = "CVC_ROOT";
    //protected static boolean initalized = false;
    protected static TreeMap<String,String> properties =null;
    protected static TreeMap<String,String> ordered = null;
    protected static HashMap<String,Element> elements = null;
    public static String ROOT_DIR = null;
    public static String CONFIG_DIR = "conf";
    public static String CONFIG_ORDER = "conf/config_order.txt";

    
    public static String getRootDir() {    	
    	if (ROOT_DIR == null) ROOT_DIR = System.getenv(ROOT_ENV);
    	if (ROOT_DIR == null) {
    		ATMLogger.error(ROOT_ENV + " not set \nEnvirontment:\n"+StaticProperties.getEnvironment());
    		throw new RuntimeException (ROOT_ENV + " not set ");
    	}
    	System.setProperty(ROOT_ENV, ROOT_DIR);
    	
    	return ROOT_DIR;
    }
    /**
     * Loads the properties found in CONFIG_FILE into the system properties.
     * Call once at start of program
     *
     * @throws JDOMException
     * @throws IOException
     */
    public static void init() throws Exception{
        getRootDir();

        StaticProperties.properties = new TreeMap<String,String>();
        StaticProperties.properties.put(StaticProperties.class.getName() + ".LOADTIME",
                new Date().toString());

        StaticProperties.elements = new HashMap<String,Element>();
        //StaticProperties.ordered = new TreeMap<String,String>();
        BufferedReader in = null;
        try {
            in = new BufferedReader(new FileReader(getConfigOrder()));
        }catch(IOException e){
            throw new Exception ("Unable to open " + getConfigOrder() ,e);
        }
        String filename= in.readLine();
        String config_orderFile = "";
        while (filename != null && filename.trim().length() != 0) {
            if (!filename.trim().startsWith("#")){
                StaticProperties.init(getConfigDir()+'/' + filename);
                config_orderFile+=getConfigDir()+'/' + filename +", ";
            }
            filename = in.readLine();
        }
        StaticProperties.properties.put(StaticProperties.class.getName() + ".CONFIG_ORDER",
                config_orderFile);
    }

    public static String getConfigDir() {
        return ROOT_DIR + "/" + CONFIG_DIR;
    }

    public static String getConfigOrder() {
        return ROOT_DIR + "/" + CONFIG_ORDER;
    }

    /**
     * Returns the value specified in config.xml.  The value could be from a Property element
     * or a Text element. Example: StaticProperties.get("mhh.atm.hla.db.HLAPeptideDB.HOST")
     *
     * @param key   This should include the package path for the class
     * @return
     */
    public static  String get(String key) {
        if (StaticProperties.properties != null){
            String prop = StaticProperties.properties.get(key);
            if (prop == null)
                ATMLogger.warn("No property found for key: " + key);
            return prop;
        }else
            throw new RuntimeException
                    ("Static Properties was not initialized. Call Tools.initialize()");    
    }

    /**
     * It uses the class name of the Object o for the namespace.  Example, in HLAPeptideDB
     * call StaticProperties.get("HOST",this), is the same as calling
     * StaticProperties.get("HOST", "mhh.atm.hla.db.HLAPeptideDB") or
     * StaticProperties.get("mhh.atm.hla.db.HLAPeptideDB.HOST")
     * @param key
     * @param o
     * @return
     */
     public static String get(String key, Object o) {
        return StaticProperties.get(o.getClass().getName() + "." + key);
    }

    public static String get(String key, Class<?> c) {
        return StaticProperties.get(c.getName() + "." + key);
    }

    /**
     *  Here, the namespace can be seperately given
     *
     * Example, in HLAPeptideDB
     * call StaticProperties.get("HOST",this), is the same as calling
     * StaticProperties.get("HOST", "mhh.atm.hla.db.HLAPeptideDB") or
     * StaticProperties.get("mhh.atm.hla.db.HLAPeptideDB.HOST")
     *
     * @param key
     * @param nameSpace
     * @return
     */
    public static String get(String key, String nameSpace) {
       return StaticProperties.get(nameSpace + "." + key);
   }


    public static String showSettings() {
        StringBuilder str = new StringBuilder(100);
        for (String key: properties.keySet()){
            str.append(key);
            str.append('=');
            str.append(properties.get(key));
            str.append('\n');
        }
        return str.toString();
    }


    public static String showSettingsHTML(boolean showEnv, boolean showSysProps) {
        StringBuilder str = new StringBuilder(100);
        str.append("<table><tr><th>Key</th><th>Value</th></tr>\n");
        if (showEnv) {
            for(String key: System.getenv().keySet()){
                str.append("<tr><td>");
                str.append(key);
                str.append("</td><td>");
                str.append(System.getenv().get(key));
                str.append("</td></tr>\n");
            }
        }
        if (showSysProps) {
            for(Enumeration<?> e = System.getProperties().propertyNames(); e.hasMoreElements();){
                String key = (String)e.nextElement();
                str.append("<tr><td>");
                str.append(key);
                str.append("</td><td>");
                str.append(System.getProperty(key));
                str.append("</td></tr>\n");
            }

        }
        for (String key: properties.keySet()){
            str.append("<tr><td>");
            str.append(key);
            str.append("</td><td>");
            str.append(properties.get(key));
            str.append("</td></tr>\n");
        }

        str.append("</table>\n");
        return str.toString();
    }
    public static String showSettingsHTML() {
        return showSettingsHTML(false, false);
    }
     /**
     *
     * @param filename
     * @throws JDOMException
     * @throws IOException
     */
    @SuppressWarnings("unchecked")
	protected static void init(String filename) throws Exception{

        SAXBuilder builder = new SAXBuilder();
        Document doc = null;
        String uriFilename = new File(filename).toURI().toString();
        try {
            doc = builder.build(uriFilename);
        } catch (JDOMException e) {
        	System.err.println("Unable to parse " + uriFilename);
            throw new Exception ("Unable to parse " + uriFilename ,e);
        } catch (IOException e) {
        	System.err.println("Unable to open " + uriFilename);
            throw new Exception ("Unable to open " + uriFilename ,e);
        }

        String /*key, value,*/ name;
        for (Element e: (List<Element>)doc.getRootElement().getChildren()) {
            name = e.getName();
            if (name.equals("SystemProperties")){
                doSystemProperty(e);
            }else if (name.equals("Package")){
                doPackage(e);
            }
        }
    
           // else System.out.println("Not an element");
        
    }

    @SuppressWarnings("unchecked")
	private static void doPackage(Element _e) {
        String name;
        for (Element e: (List<Element>)_e.getChildren()){
            name = e.getName();
            if (name.equals("Package")){
                doPackage(e);
            }else doClass(e);
        }
    }

    @SuppressWarnings("unchecked")
	private static void doClass(Element _e) {
        String prefix = _e.getNamespaceURI() + '.' + _e.getName();
        String childName;
        Attribute type;


        for (Element e: (List<Element>)_e.getChildren()){
            childName = e.getName();
            if (childName.equals("Properties")){
                for (Element property : (List<Element>)e.getChildren()){
                    type = property.getAttribute("type");
                    if (type != null && type.getValue().equals("relativePath")){
                        StaticProperties.properties.put(prefix + "." +
                            property.getName(), ROOT_DIR + "/" + property.getText());

                    }else{
                        StaticProperties.properties.put(prefix + "." +
                            property.getName(), property.getText());

                    }
                }
            }else if (childName.equals("ElementProperties")){
                 for (Element property: (List<Element>)e.getChildren()){
                     StaticProperties.elements.put(prefix + "." +
                            property.getName(), property);
                     
                }
            }
             System.setProperty(e.getName(), e.getText());
        }
    }

    @SuppressWarnings("unchecked")
	private static void doSystemProperty(Element _e) {
        for (Element e: (List<Element>)_e.getChildren()){
             System.setProperty(e.getName(), e.getText());
        }
    }


    public static void main (String[] args) {
        try {
            Tools.initialize();

            System.out.println("Config folder: " + StaticProperties.getConfigDir());
            System.out.println("ENVIRONMENT");
            System.out.println(StaticProperties.getEnvironment());

            System.out.println("SYSTEM PROPERTIES");
            System.out.println(StaticProperties.getSystemProperties());
         } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    public static String getEnvironment() {
        StringBuilder str = new StringBuilder(100);
        str.append("ENVIRONMENT\n");
        for(String key: System.getenv().keySet()){
            str.append(key);
            str.append (" = ");
            str.append(System.getenv().get(key));
            str.append('\n');
        }
        return str.toString();
    }

    public static String getSystemProperties() {
        StringBuilder str = new StringBuilder(100);
        str.append("PROPERTIES\n");
        for(Enumeration<?> e = System.getProperties().propertyNames(); e.hasMoreElements();){
            String key = (String)e.nextElement();
            str.append(key);
            str.append (" = ");
            str.append(System.getProperty(key));
            str.append('\n');
        }
        return str.toString();

    }

    public static boolean isMultiProcessor() {
        return !System.getenv("NUMBER_OF_PROCESSORS").equals("1");
    }
    
    
}

