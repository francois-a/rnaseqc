package org.broadinstitute.cga.rnaseq;

/**
 * Created by IntelliJ IDEA.
 * User: ddeluca
 * Date: 3/25/11
 * Time: 2:05 PM
 */
public class Annotations {

    String type;
    private String strand;
    String geneId;

    private String index ="";
    private String sequence ="";
    private String name;
    private Float gCContent = null;

    @Override
    public String toString() {
        return geneId +" " + name + " " + type + " " + strand;
    }

    public void setType(String typ) {
        this.type = typ;
    }

    public boolean isReverseStrand() {
        return this.strand.equals("-1") || this.strand.equals("-");
    }

    public void setGeneId(String id) {
        this.geneId = id;
    }

    public String getGeneId() {
        return this.geneId;
    }

    public String getName(){
        return name;
    }

    public void setIndex(String index) {
        this.index = index;
    }

    public void setSequence(String seq) {
        this.sequence = seq;
    }

    public String getSequence() {
        return sequence;
    }

    public float getGCContent() {
        if (this.gCContent == null) {
            this.gCContent = getGCContent(this.sequence);
        }
        return this.gCContent;
    }

    public static float getGCContent(String sequence) {
        int gcCount = 0;
        for (char c: sequence.toCharArray()) {
            if (c == 'G' || c=='g' || c=='C' || c=='c') {
                gcCount++;
            }
        }
        return (float)gcCount / (float)sequence.length();
    }

    public String getType() {
        return type;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public String getIndex() {
        return index;
    }

    public String getStrand() {
        return strand;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setGCContent(Float GCContent) {
        this.gCContent = GCContent;
    }
}
