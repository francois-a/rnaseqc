package org.broadinstitute.cga.tools;

import org.apache.log4j.Logger;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.xml.DOMConfigurator;


/**
 * User: DeLuca
 * Date: 19.01.2006
 * Time: 12:08:09
 */
public class ATMLogger {

    public static void init(String configFileName, String rootPath) {

        ATMLogger.logger = Logger.getLogger(ATMLogger.class);
        if (configFileName == null) {
            BasicConfigurator.configure();
        } else {
            try {
                DOMConfigurator.configure(configFileName);
            } catch (Exception e) {
                throw new RuntimeException("Could not configure Logger for file " +
                        configFileName + ". Reason: " + e.getMessage());
            }
        }
        StaticProperties.ROOT_DIR = rootPath;
        ATMLogger.logger.debug("log4j config file: " + configFileName);
        ATMLogger.logger.debug("logFiles file: " + configFileName);
    }

    public static void changeLogger(String loggerName) {
        ATMLogger.logger = Logger.getLogger(loggerName);
    }

    public static Logger getLogger(String name) {
        return Logger.getLogger(name);
    }

    public static Logger getLogger(){
        return ATMLogger.logger;
    }

    public static void debug(Object message){
        if (logger == null)
            System.out.println(message);
        else ATMLogger.logger.debug(message);
    }
    
    public static void info(Object message){
        if (logger == null)
            System.out.println(message);
        else ATMLogger.logger.info(message);
    }
    
    public static void warn(Object message){
        if (logger == null)
            System.out.println(message);
        else ATMLogger.logger.warn(message);
    }
    
    public static void error(Object message){
        if (logger == null)
            System.out.println(message);
        else ATMLogger.logger.error(message);
    }
    
    public static void fatal(Object message){
        if (logger == null)
            System.out.println(message);
        else ATMLogger.logger.fatal(message);
    }

    public static void debug(Object message,  Throwable t){
        if (logger == null)
            System.out.println(message);
        else ATMLogger.logger.debug( message, t);
    }

    public static void info(Object message,  Throwable t){
        if (logger == null)
            System.out.println(message);
        else ATMLogger.logger.info( message, t);
    }

    public static void warn(Object message,  Throwable t){
        if (logger == null)
            System.out.println(message);
        else ATMLogger.logger.warn( message, t);
    }

    public static void error(Object message,  Throwable t){
        if (logger == null)
            System.out.println(message);
        else ATMLogger.logger.error( message, t);
    }

    public static void fatal(Object message,  Throwable t){
        if (logger == null)
            System.out.println(message);
        else ATMLogger.logger.fatal( message, t);
    }

    protected static Logger logger = null;

}
