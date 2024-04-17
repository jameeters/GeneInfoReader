package org.pankratzlab;

import java.util.HashSet;
import java.util.Set;

public class GeneGrouping {
  private final Set<BasicFeature> genes = new HashSet<>();
  private final Set<BasicFeature> mainContigGenes = new HashSet<>();
  public final String geneId;

  public GeneGrouping(String geneId) {
    this.geneId = geneId;
  }

  public GeneGrouping(String geneId, Iterable<BasicFeature> genes){
    this(geneId);
    for(BasicFeature gene : genes){
      this.addGene(gene);
    }
  }

  public void addGene(BasicFeature gene) {
    if (!gene.isGene()) {
      throw new IllegalArgumentException("This is not a gene, you can't add it to a gene group");
    }
    if (!gene.xRefGeneId.equals(this.geneId)) {
      throw new IllegalArgumentException(
          "This gene (" + gene.id + ") has the wrong xRefGeneId: " + gene.xRefGeneId + " expected "
          + this.geneId);
    }
    this.genes.add(gene);
    if (gene.onMainContig) {
      this.mainContigGenes.add(gene);
    }
  }

  public boolean hasMainContigGene() {
    return this.mainContigGenes.size() > 0;
  }

  public boolean hasMultipleMainContigGenes() {
    return this.mainContigGenes.size() > 1;
  }

  public Set<BasicFeature> getGenes() {
    return this.genes;
  }

  public Set<BasicFeature> getMainContigGenes(){
    return this.mainContigGenes;
  }

  public int countTotalGenes() {
    return this.genes.size();
  }

  public int countMainContigGenes() {
    return this.mainContigGenes.size();
  }

  public int getChr() {
    if (this.hasMainContigGene()) {
      return mainContigGenes.stream().findAny().get().getChr();
    } else {
      return genes.stream().findAny().get().getChr();
    }
  }

  public int getStart() {
    if (this.hasMainContigGene()) {
      return mainContigGenes.stream().findAny().get().start;
    } else {
      return genes.stream().findAny().get().start;
    }
  }

  public int getEnd() {
    if (this.hasMainContigGene()) {
      return mainContigGenes.stream().findAny().get().end;
    } else {
      return genes.stream().findAny().get().end;
    }
  }

  public int compareTo (GeneGrouping other) {
    if (other.getChr() != this.getChr()) {
      return Integer.compare(this.getChr(), other.getChr());
    }
    if (other.getStart() != this.getStart()) {
      return Integer.compare(this.getStart(), other.getStart());
    }
    return Integer.compare(this.getEnd(), other.getEnd());
  }
}
