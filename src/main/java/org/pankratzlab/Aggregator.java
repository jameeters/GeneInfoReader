package org.pankratzlab;

import htsjdk.tribble.gff.Gff3BaseData;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class Aggregator {
  final Map<String, BasicFeature> featureMap = new HashMap<>();
  final Set<BasicFeature> genes = new HashSet<>();

  public Aggregator() {
    //no op
  }

  public void add(Gff3BaseData baseData) {
    BasicFeature feat = new BasicFeature(baseData);
    this.featureMap.put(feat.id, feat);
    try {
      String parentId = baseData.getAttribute("Parent").get(0);
      feat.parent = this.featureMap.get(parentId);
      feat.parent.children.add(feat);
    } catch (IndexOutOfBoundsException e) {
      feat.parent = null;
    }
  }

  void findGenesAndExons() {
    for (BasicFeature feat : featureMap.values()) {
      if (feat.type.equals("gene")) {
        this.genes.add(feat);
        feat.getDescendantExons();
      }
    }
  }

  void writeTestOutput() throws IOException {
    String fileName = "/tmp/geneinfo";
    FileWriter writer = new FileWriter(fileName);

    for (BasicFeature gene : genes) {
      writer.write(gene.id + "  " + gene.start + "  " + gene.end + "\n");
      for (BasicFeature exon : gene.getDescendantExons()) {
        writer.write("   |" + exon.id + "  " + exon.start + "  " + exon.end + "\n");
      }
    }
    writer.close();
  }

  void writeQcOutput() throws IOException {
    String fileName = "/tmp/geneinfoQC";
    FileWriter writer = new FileWriter(fileName);

    int[] chrGeneCounts = new int[25];

    for (BasicFeature gene : genes) {
      chrGeneCounts[gene.getChr()]++;
    }
    writer.write("chr\tgeneCount\n");
    for (int i = 0; i < chrGeneCounts.length; i++){
      writer.write(i + "\t" + chrGeneCounts[i] + "\n");
    }
      writer.close();
  }
}
