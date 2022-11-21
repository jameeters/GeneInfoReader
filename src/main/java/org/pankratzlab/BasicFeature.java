package org.pankratzlab;

import htsjdk.tribble.gff.Gff3BaseData;

import java.util.HashSet;
import java.util.Set;

public class BasicFeature {
  private static final Aggregator aggregator = new Aggregator();

  final String id, type;
  private BasicFeature parent;
  final Set<BasicFeature> children = new HashSet<>();

  boolean exonsFound = false;
  Set<BasicFeature> descendantExons = new HashSet<>();

  public BasicFeature(Gff3BaseData baseData) {
    this.type = baseData.getType();
    this.id = baseData.getId();
    aggregator.featureMap.put(this.id, this);
    try {
      String parentId = baseData.getAttribute("Parent").get(0);
      this.parent = aggregator.featureMap.get(parentId);
      this.parent.children.add(this);
    } catch (IndexOutOfBoundsException e) {
      this.parent = null;
    }

  }

  public Set<BasicFeature> getDescendantExons() {
    if(this.exonsFound) {
      return this.descendantExons;
    }
    Set<BasicFeature> exons = new HashSet<>();
    for (BasicFeature child : this.children) {
      exons.addAll(child.getDescendantExons());
    }

    if(this.type.equals("exon") && exons.size() > 0) {
      System.out.println("hey this exon has descendant exons, isn't that weird?!");
    }

    if(this.type.equals("exon")) {
      exons.add(this);
    }

    this.descendantExons = exons;
    this.exonsFound = true;
    return exons;
  }

  static int count() {
    return aggregator.featureMap.size();
  }
}
