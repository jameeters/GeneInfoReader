package org.pankratzlab;

import htsjdk.tribble.gff.Gff3BaseData;
import org.pankratzlab.common.filesys.GeneData;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class BasicFeature {

  final String id, type;
  BasicFeature parent;
  final Set<BasicFeature> children = new HashSet<>();

  boolean exonsFound = false;
  Set<BasicFeature> descendantExons = new HashSet<>();

  final int start, end;
  final String name;
  final String contig;
  final byte strand;

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
        // todo: is this right?
        tempStrand = GeneData.BOTH_STRANDS;
    }
    this.strand = tempStrand;
  }

  public Set<BasicFeature> getDescendantExons() {
    if (this.exonsFound) {
      return this.descendantExons;
    }
    Set<BasicFeature> exons = new HashSet<>();
    for (BasicFeature child : this.children) {
      exons.addAll(child.getDescendantExons());
    }

    if (this.type.equals("exon") && exons.size() > 0) {
      System.out.println("hey this exon has descendant exons, isn't that weird?!");
    }

    if (this.type.equals("exon")) {
      exons.add(this);
    }

    this.descendantExons = exons;
    this.exonsFound = true;
    return exons;
  }

  public byte getChr() {
    switch (contig) {
      case "NC_000001.11":
        return 1;
      case "NC_000002.12":
        return 2;
      case "NC_000003.12":
        return 3;
      case "NC_000004.12":
        return 4;
      case "NC_000005.10":
        return 5;
      case "NC_000006.12":
        return 6;
      case "NC_000007.14":
        return 7;
      case "NC_000008.11":
        return 8;
      case "NC_000009.12":
        return 9;
      case "NC_000010.11":
        return 10;
      case "NC_000011.10":
        return 11;
      case "NC_000012.12":
        return 12;
      case "NC_000013.11":
        return 13;
      case "NC_000014.9":
        return 14;
      case "NC_000015.10":
        return 15;
      case "NC_000016.10":
        return 16;
      case "NC_000017.11":
        return 17;
      case "NC_000018.10":
        return 18;
      case "NC_000019.10":
        return 19;
      case "NC_000020.11":
        return 20;
      case "NC_000021.9":
        return 21;
      case "NC_000022.11":
        return 22;
      case "NC_000023.11":
        return 23;
      case "NC_000024.10":
        return 24;
      default:
        return 0;

    }
  }

  public int[] getBoundariesAsArray() {
    return new int[] {start, end};
  }

  public int[][] getDescendantExonBoundariesAsArray() {
    //todo: memoize
    this.getDescendantExons();
    int[][] bounds = new int[descendantExons.size()][];
    List<BasicFeature> exonsList = new ArrayList<>(descendantExons);
    for (int i = 0; i < exonsList.size(); i++) {
      bounds[i] = exonsList.get(i).getBoundariesAsArray();
    }
    return bounds;
  }

  public GeneData toGeneData() {
    //todo memoize
    if (!this.type.equals("gene")) {
      throw new IllegalStateException("GeneData should only be created onn genes");
    }

    String[] ncbiAssessionNums = new String[0];

    //todo: positionFinalized, multiLoc, collapsedIsoFormGene still need to be done correctly
    GeneData geneData = new GeneData(name, ncbiAssessionNums, this.getChr(), true, strand, start,
                                     end, this.getDescendantExonBoundariesAsArray(), (byte) 0,
                                     false);
    return geneData;
  }

}
