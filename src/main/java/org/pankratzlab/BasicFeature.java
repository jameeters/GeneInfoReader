package org.pankratzlab;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import org.pankratzlab.common.filesys.GeneData;

import htsjdk.tribble.gff.Gff3BaseData;
import htsjdk.tribble.gff.Gff3Feature;

public class BasicFeature {

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

  private BasicFeature(String id, String type, int start, int end, String name, String contig,
                       byte strand, String xRefGeneId) {
    this.type = type;
    this.id = id;
    this.start = start;
    this.end = end;
    this.name = name;
    this.contig = contig;
    this.onMainContig = contigToChrMapping.containsKey(this.contig);
    this.strand = strand;
    this.xRefGeneId = xRefGeneId;
  }

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
    if (!this.exonsFound) {
      this.getDescendantExons();
    }
    Set<BasicFeature> introns = new HashSet<>();

    List<BasicFeature> orderedExons = this.descendantExons.stream().filter(distinctByLocation()).sorted(BasicFeature::compareLocation)
        .collect(Collectors.toList());

    // Introns are everything in the gene that's not an exon
    // start by looking between the start of the gene and start of the first exon
    int prevEnd = this.start;
    for (BasicFeature exon : orderedExons) {
      if (exon.start > prevEnd) {
        // new thing
        String id = this.id + "_intron_" + start + "_" + end;
        // we've already made sure that prevEnd is outside an exon
        int intronStart = prevEnd;
        // the intron and exon don't share a base pair
        int intronEnd = exon.start - 1;
        BasicFeature intron = new BasicFeature(id, "intron", intronStart, intronEnd, id, this.contig, this.strand, this.xRefGeneId);
        introns.add(intron);
      }
      prevEnd = exon.end + 1;
    }
    if (this.end > prevEnd) {
      // one more intron between the last exon and the end of the gene
      BasicFeature intron = new BasicFeature(id, "intron", prevEnd, this.end, id, this.contig, this.strand, this.xRefGeneId);
      introns.add(intron);
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

  public String getPositionAsString() {
    return String.join("_", String.valueOf(getChr()), String.valueOf(start), String.valueOf(end));
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

  public int compareLocation (BasicFeature other) {
    if (other.getChr() != this.getChr()) {
      return Integer.compare(this.getChr(), other.getChr());
    }
    if (other.start != this.start) {
      return Integer.compare(this.start, other.start);
    }
    return Integer.compare(this.end, other.end);
  }

  public List<String> toBedLines(boolean includeExons, boolean includeIntrons) {
    if (!this.isGene()) {
      throw new RuntimeException("I'm not a gene and I don't want to be turned into bed lines!");
    }
    if (getChr() == 26) {
      return List.of();
    }
    List<BasicFeature> children = new ArrayList<>();
    if (includeExons) {
      children.addAll(this.descendantExons.stream().filter(distinctByLocation()).collect(Collectors.toList()));
    }
    if (includeIntrons) {
      children.addAll(this.descendantIntrons.stream().filter(distinctByLocation()).distinct().collect(Collectors.toList()));
    }

    children.sort(BasicFeature::compareLocation);
    List<String> lines =  new ArrayList<>();
    for (int i = 0; i < children.size(); i++) {
      BasicFeature child = children.get(i);

      String name = String.join("_", this.name, child.type.substring(0,1), String.valueOf(i));

      String chrstr = String.valueOf(child.getChr());
      if (child.getChr() == 23) chrstr = "X";
      if (child.getChr() == 24) chrstr = "Y";

      lines.add(String.join("\t",
                            "chr" + chrstr,
                            String.valueOf(child.start),
                            String.valueOf(child.end),
                            name));
    }
    return lines;
  }

  public static Predicate<BasicFeature> distinctByLocation() {
    Set<Object> seen = ConcurrentHashMap.newKeySet();
    return t -> seen.add(t.getPositionAsString());
  }
}
