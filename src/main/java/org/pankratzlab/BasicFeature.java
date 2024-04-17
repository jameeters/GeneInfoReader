package org.pankratzlab;

import htsjdk.tribble.gff.Gff3BaseData;
import htsjdk.tribble.gff.Gff3Feature;
import org.pankratzlab.common.filesys.GeneData;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class BasicFeature implements Comparable<BasicFeature> {

  private static final Map<String, Integer> contigToChrMapping = Map.ofEntries(
      new AbstractMap.SimpleEntry<>("NC_000001.11", 1),
      new AbstractMap.SimpleEntry<>("NC_000002.12", 2),
      new AbstractMap.SimpleEntry<>("NC_000003.12", 3),
      new AbstractMap.SimpleEntry<>("NC_000004.12", 4),
      new AbstractMap.SimpleEntry<>("NC_000005.10", 5),
      new AbstractMap.SimpleEntry<>("NC_000006.12", 6),
      new AbstractMap.SimpleEntry<>("NC_000007.14", 7),
      new AbstractMap.SimpleEntry<>("NC_000008.11", 8),
      new AbstractMap.SimpleEntry<>("NC_000009.12", 9),
      new AbstractMap.SimpleEntry<>("NC_000010.11", 10),
      new AbstractMap.SimpleEntry<>("NC_000011.10", 11),
      new AbstractMap.SimpleEntry<>("NC_000012.12", 12),
      new AbstractMap.SimpleEntry<>("NC_000013.11", 13),
      new AbstractMap.SimpleEntry<>("NC_000014.9", 14),
      new AbstractMap.SimpleEntry<>("NC_000015.10", 15),
      new AbstractMap.SimpleEntry<>("NC_000016.10", 16),
      new AbstractMap.SimpleEntry<>("NC_000017.11", 17),
      new AbstractMap.SimpleEntry<>("NC_000018.10", 18),
      new AbstractMap.SimpleEntry<>("NC_000019.10", 19),
      new AbstractMap.SimpleEntry<>("NC_000020.11", 20),
      new AbstractMap.SimpleEntry<>("NC_000021.9", 21),
      new AbstractMap.SimpleEntry<>("NC_000022.11", 22),
      new AbstractMap.SimpleEntry<>("NC_000023.11", 23),
      new AbstractMap.SimpleEntry<>("NC_000024.10", 24),
      new AbstractMap.SimpleEntry<>("NC_012920.1", 26));

  final String id, type;
  BasicFeature parent;
  final Set<BasicFeature> children = new HashSet<>();

  boolean exonsFound = false;
  Set<BasicFeature> descendantExons = new HashSet<>();

  boolean intronsFound = false;
  Set<BasicFeature> descendantIntrons = new HashSet<>();

  final int start, end;
  final String name;
  final String contig;
  final byte strand;

  final boolean onMainContig;
  final String xRefGeneId;

  private int[][] descendantExonBoundaries;
  private GeneData geneData;

  public BasicFeature(Gff3BaseData baseData) {
    this.type = baseData.getType();
    this.id = baseData.getId();
    this.start = baseData.getStart();
    this.end = baseData.getEnd();
    this.name = baseData.getName();
    this.contig = baseData.getContig();

    byte tempStrand = -1;
    switch (baseData.getStrand().encodeAsChar()) {
      case '+':
        tempStrand = GeneData.PLUS_STRAND;
        break;
      case '-':
        tempStrand = GeneData.MINUS_STRAND;
        break;
      case '.':
        tempStrand = GeneData.BOTH_STRANDS;
    }
    this.strand = tempStrand;
    this.onMainContig = contigToChrMapping.containsKey(this.contig);
    this.xRefGeneId = findXRefGeneId(baseData);
  }

  public BasicFeature(Gff3Feature feature) {
    this(feature.getBaseData());
  }

  private static String findXRefGeneId(Gff3BaseData baseData) {
    List<String> provisionalDbXRef = baseData.getAttribute("Dbxref");
    for (String ref : provisionalDbXRef) {
      if (ref.startsWith("GeneID")) {
        return ref;
      }
    }
    return Aggregator.BAD_OR_MISSING;
  }

  public Set<BasicFeature> getDescendantExons() {
    if (this.exonsFound) {
      return this.descendantExons;
    }
    Set<BasicFeature> exons = new HashSet<>();
    for (BasicFeature child : this.children) {
      exons.addAll(child.getDescendantExons());
    }

    if (this.isExon() && exons.size() > 0) {
      System.out.println("hey this exon has descendant exons, isn't that weird?!");
    }

    if (this.isExon()) {
      exons.add(this);
    }

    this.descendantExons = exons;
    this.exonsFound = true;
    return exons;
  }

  public Set<BasicFeature> getDescendantIntrons() {
    if (this.intronsFound) {
      return this.descendantIntrons;
    }
    Set<BasicFeature> introns = new HashSet<>();
    for (BasicFeature child : this.children) {
      introns.addAll(child.getDescendantIntrons());
    }

    if (this.isIntron() && introns.size() > 0) {
      System.out.println("hey this intron has descendant introns, isn't that weird?!");
    }

    if (this.isIntron()) {
      introns.add(this);
    }

    this.descendantIntrons = introns;
    this.intronsFound = true;
    return introns;
  }

  public byte getChr() {
    return contigToChrMapping.getOrDefault(this.contig, 0).byteValue();
  }

  public int[] getBoundariesAsArray() {
    return new int[] {start, end};
  }

  public int[][] getDescendantExonBoundariesAsArray() {
    if (descendantExonBoundaries == null) {
      this.getDescendantExons();
      descendantExonBoundaries = new int[descendantExons.size()][];
      List<BasicFeature> exonsList = new ArrayList<>(descendantExons);
      for (int i = 0; i < exonsList.size(); i++) {
        descendantExonBoundaries[i] = exonsList.get(i).getBoundariesAsArray();
      }
    }
    return descendantExonBoundaries;
  }

  public GeneData toGeneData() {
    if (geneData == null) {
      if (!this.isGene()) {
        throw new IllegalStateException("GeneData should only be created on genes");
      }
      // todo: ncbi numbers
      String[] ncbiAssessionNums = new String[0];
      //todo: positionFinalized, multiLoc, collapsedIsoFormGene still need to be done correctly
      geneData = new GeneData(name, ncbiAssessionNums, this.getChr(), true, strand, start, end,
                              this.getDescendantExonBoundariesAsArray(), (byte) 0, false);
    }
    return geneData;
  }

  public boolean isGene(){
    return this.type.equals("gene");
  }

  public boolean isExon() {
    return this.type.equals("exon");
  }

  public boolean isIntron() {
    return this.type.equals("intron");
  }

  public int compareTo (BasicFeature other) {
    if ((other.parent != null && this.parent != null) && other.parent != this.parent) {
      return this.parent.compareTo(other.parent);
    }
    if (other.getChr() != this.getChr()) {
      return this.getChr() - other.getChr();
    }
    if (other.start != this.start) {
      return this.start - other.start;
    }
    return this.end - other.end;
  }

  public List<String> toBedLines() {
    if (!this.isGene()) {
      throw new RuntimeException("I'm not a gene and I don't want to be turned into bed lines!");
    }
    List<BasicFeature> children = new ArrayList<>();
    if (exonsFound) {
      children.addAll(this.descendantExons);
    }
    if (intronsFound) {
      children.addAll(this.descendantIntrons);
    }
    children.sort(BasicFeature::compareTo);
    List<String> lines =  new ArrayList<>();
    for (int i = 0; i < children.size(); i++) {
      BasicFeature child = children.get(i);
      String name = String.join("_", this.name, child.type.substring(0,1), String.valueOf(i));
      lines.add(String.join("\t",
                            String.valueOf(child.getChr()),
                            String.valueOf(child.start),
                            String.valueOf(child.end),
                            name));
    }
    return lines;
  }
}
