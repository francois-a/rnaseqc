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

 */

public class Tools {

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


    /**
     * Performs the initializations necessary when an application starts.
     * 1) Initializes the StaticProperties class using the conf/config_order.txt file
     * 2) Initializes the ATMLogger using log4j xml file, as defined in static properties
     */
    static public void initialize(){
        try {
            if (StaticProperties.properties == null) {
                StaticProperties.init();
                ATMLogger.init(StaticProperties.get("LOG4J_CONF_FILE", ATMLogger.class),
                    StaticProperties.ROOT_DIR );
                ATMLogger.debug(StaticProperties.showSettings());
                ATMLogger.debug("ATM_ROOT: " + StaticProperties.ROOT_DIR);
            }
        } catch (Exception e) {
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


    public static String stripFileName(String filePath) {
        File f = new File(filePath);
        return f.getName();
    }

}
