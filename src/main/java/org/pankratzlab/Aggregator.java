package org.pankratzlab;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class Aggregator {
  final Map<String, BasicFeature> featureMap = new HashMap<>();
  final Set<BasicFeature> genes = new HashSet<>();

  public Aggregator(){
    //no op
  }

  private void findGenes() {
    for (BasicFeature feat : featureMap.values()) {
      if(feat.type.equals("gene")){
        this.genes.add(feat);
        feat.getDescendantExons();
      }
    }
  }


}
