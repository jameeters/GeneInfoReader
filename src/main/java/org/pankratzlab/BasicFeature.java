package org.pankratzlab;

import htsjdk.tribble.gff.Gff3BaseData;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class BasicFeature {
  private static final Map<String, BasicFeature> featureMap = new HashMap<>();
  private static final Map<String, BasicFeature> geneMap = new HashMap<>();

  final String id, type;
  private BasicFeature parent;
  final Set<BasicFeature> children = new HashSet<>();

  public BasicFeature(Gff3BaseData baseData) {
    this.type = baseData.getType();
    this.id = baseData.getId();
    featureMap.put(this.id, this);
    try {
      String parentId = baseData.getAttribute("Parent").get(0);
      this.parent = featureMap.get(parentId);
      this.parent.children.add(this);
    } catch (IndexOutOfBoundsException e) {
      this.parent = null;
    }

  }

  static int count() {
    return featureMap.size();
  }
}
