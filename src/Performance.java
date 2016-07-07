package org.broadinstitute.cga.tools;

/**
 * This class works as a timer for measuring performance
 * @author David
 *
 */
public class Performance {

    public enum Resolution {nanoseconds, milliseconds, seconds, minutes};
    protected Resolution res;
    long started;
    String title;
    public Performance() {
        this(null,Resolution.seconds);
    }

    public Performance(Resolution resolution) {
        this(null,resolution);
    }

    /**
     * The default resolution is seconds
     * @param name
     */
    public Performance(String name) {
        this(name,Resolution.seconds);		
    }

    public Performance(String name, Resolution resolution) {
        this.res = resolution;
        this.title = name;
        this.started = System.nanoTime();
    }

    @Override
    public String toString() {
        return title + ":\t" + this.getTimeDiffString();
    }

    public void start() {
        started=System.nanoTime();
    }

    public long getExtactNanosTook() {
        return System.nanoTime()-started;
    }

    /**
     * Gets the time took in seconds
     * @return
     */
	public long getTimeDiff() {
        switch (this.res) {
            case nanoseconds:
                return  System.nanoTime()-started ;
            case milliseconds:
                return (System.nanoTime()-started) / 1000000l;
            case seconds:
                return (System.nanoTime()-started) / 1000000000l;
            case minutes:
                return (System.nanoTime()-started) / 60000000000l;
        }
        throw new RuntimeException ("Resolution: " + this.res + " not found in enum ");
    }

    /**
     * Gets the time took in seconds
     * @return
     */
    public String getTimeDiffString() {
        switch (this.res) {
            case nanoseconds:
                return  String.valueOf(System.nanoTime()-started) +" ns" ;
            case milliseconds:
                return  String.valueOf( (System.nanoTime()-started) / 1000000l)+" ms";
            case seconds:
                return String.valueOf(((System.nanoTime()-started) / 1000000000l))+" s";
            case minutes:
                return String.valueOf((System.nanoTime()-started) / 60000000000l) +" min" ;
        }
        throw new RuntimeException ("Resolution: " + this.res + " not found in enum ");
    }

    public void outputIfTookTime() {
        if (this.getTimeDiff() > 0) {
            System.out.println(this.toString());
        }
    }

    public static void main (String[]args){
        Performance p = new Performance("my perf", Resolution.milliseconds);
        for (int i=0; i < 10000;i++) {}
        System.out.println(p);
    }
    
}
