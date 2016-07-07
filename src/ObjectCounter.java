package org.broadinstitute.cga.tools;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Set;
import java.util.StringTokenizer;

public class ObjectCounter  <T> {

    int total = 0;

    HashMap<T,Integer> objToInt = new HashMap<T,Integer>();

    public ObjectCounter(){}

    public ObjectCounter(String str) {
        StringTokenizer toks = new StringTokenizer(str, "(),");

        while (toks.hasMoreTokens()) {
            T obj = (T)toks.nextElement();
            Integer i = Integer.valueOf(toks.nextToken());
            total += i;
            objToInt.put(obj, i);
        }
    }


    public void add(T obj) {
        total ++;
        Integer i = objToInt.get(obj);
        if (i == null) {
            i = new Integer(1);
            objToInt.put(obj, i);
        } else {
            objToInt.put(obj, new Integer( i+1));
        }
    }


    public int getTotal() {
        return total;
    }


    @Override
    public String toString() {
        StringBuilder str = new StringBuilder();
        if (total == 0) {
            return "";
        }
        int c = 0;
        Set<T> set = objToInt.keySet();
        int size = set.size();
        for (T key: set) {
            c++;
            str.append(key).append("(").append(objToInt.get(key));
            if (c == size) {
                str.append(")");
            } else {
                str.append("),");
            }
        }
        return str.toString();
    }


    public int getNumCounters() {
        return objToInt.size();
    }


    public Set<T> getKeySet() {
        return objToInt.keySet();
    }


    public int getCount(T key) {
        Integer count = objToInt.get(key);
        if (count == null) {
            return 0;
        } else {
            return count;
        }
    }


    public void printTable(PrintStream out) {

        if (total == 0) {
            out.println( "Empty");
        }
        int c = 0;
        Set<T> set = objToInt.keySet();
        for (T key: set) {
            c++;
            out.print(key+"\t"+objToInt.get(key)+"\n");
        }
    }

}
