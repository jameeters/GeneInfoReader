package org.pankratzlab;

import htsjdk.tribble.gff.Gff3BaseData;
import htsjdk.tribble.gff.Gff3Feature;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

public class Aggregator {
  final Map<String, BasicFeature> featureMap = new HashMap<>();
  final Set<BasicFeature> genes = new HashSet<>();
  final Set<String> duplicateIds = new HashSet<>();

  final static String BAD_OR_MISSING = "BAD_OR_MISSING";
  Map<String, Set<BasicFeature>> genesByXRefGeneId = new TreeMap<>();

  public Aggregator() {
    //no op
  }

  private void add(BasicFeature feat, Gff3BaseData baseData) {

    if (featureMap.containsKey(feat.id)) {
      duplicateIds.add(feat.id);
    }
    this.featureMap.put(feat.id, feat);
    try {
      String parentId = baseData.getAttribute("Parent").get(0);
      feat.parent = this.featureMap.get(parentId);
      feat.parent.children.add(feat);
    } catch (IndexOutOfBoundsException e) {
      feat.parent = null;
    }
  }

  public void add(Gff3BaseData baseData) {
    BasicFeature feat = new BasicFeature(baseData);
    this.add(feat, baseData);
  }

  public void add(Gff3Feature superFeature) {
    BasicFeature feat = new BasicFeature(superFeature);
    this.add(feat, superFeature.getBaseData());
  }

  void findGenesAndExons() {
    for (BasicFeature feat : featureMap.values()) {
      if (feat.type.equals("gene")) {
        this.genes.add(feat);
        feat.getDescendantExons();
      }
    }
  }

  public void computeXRefMap() {
    for (BasicFeature gene : genes) {
      genesByXRefGeneId.computeIfAbsent(gene.xRefGeneId, (String xRefId) -> genes.parallelStream().filter(
          g -> g.xRefGeneId.equals(xRefId)).collect(Collectors.toSet()));
      if(genesByXRefGeneId.size()%2000 == 0){
        System.out.println(genesByXRefGeneId.size() + " entries computed");
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
    String chrGeneCountsFile = "/tmp/chrGeneCounts.tsv";
    String seqIdCountsFile = "/tmp/seqIdCounts.tsv";
    String duplicateIdsFile = "/tmp/duplicateIds.tsv";
    String genesContigsFile = "/tmp/genesContigs.tsv";
    String geneIdMappingFile = "/tmp/geneIdMapping.tsv";

    List<String> fileNames = List.of(chrGeneCountsFile, seqIdCountsFile, duplicateIdsFile,
                                     genesContigsFile, geneIdMappingFile);

    for (String fileName : fileNames) {
      File f = new File(fileName);
      boolean result = Files.deleteIfExists(f.toPath());
      if (result) {
        System.out.println("deleted " + fileName);
      }
    }

    FileWriter geneCountsWriter = new FileWriter(chrGeneCountsFile);
    FileWriter seqIdCountsWriter = new FileWriter(seqIdCountsFile);

    int[] chrGeneCounts = new int[27];
    Map<String, Integer> seqIdCounts = new TreeMap<>();
    Map<String, Integer> seqIdTochrMapping = new TreeMap<>();

    for (BasicFeature gene : genes) {
      chrGeneCounts[gene.getChr()]++;
      if (seqIdCounts.containsKey(gene.contig)) {
        Integer newCount = seqIdCounts.get(gene.contig) + 1;
        seqIdCounts.put(gene.contig, newCount);
      } else {
        seqIdCounts.put(gene.contig, 1);
        seqIdTochrMapping.put(gene.contig, (int) gene.getChr());
      }
    }
    geneCountsWriter.write("chr\tgeneCount\n");
    for (int i = 0; i < chrGeneCounts.length; i++) {
      geneCountsWriter.write(i + "\t" + chrGeneCounts[i] + "\n");
    }
    geneCountsWriter.close();

    seqIdCountsWriter.write("seqId\tgeneCount\tchrMapping\n");
    for (Map.Entry<String, Integer> entry : seqIdCounts.entrySet()) {
      String contig = entry.getKey();
      Integer count = entry.getValue();
      Integer chrMapping = seqIdTochrMapping.get(contig);
      seqIdCountsWriter.write(contig + "\t" + count + "\t" + chrMapping + "\n");
    }
    seqIdCountsWriter.close();

    FileWriter duplicateIdsWriter = new FileWriter(duplicateIdsFile);
    duplicateIdsWriter.write("id\tinGenes\tinExons\n");

    Set<String> geneIds = genes.stream().map(gene -> gene.id).collect(Collectors.toSet());
    Set<String> exonIds = genes.stream().flatMap(gene -> gene.descendantExons.stream())
                               .map(exon -> exon.id).collect(Collectors.toSet());

    for (String id : duplicateIds) {
      int duplicateGene = geneIds.contains(id) ? 1 : 0;
      int duplicateExon = exonIds.contains(id) ? 1 : 0;
      duplicateIdsWriter.write(id + "\t" + duplicateGene + "\t" + duplicateExon + "\n");
    }
    duplicateIdsWriter.close();

    FileWriter geneContigWriter = new FileWriter(genesContigsFile);
    geneContigWriter.write("id\tcontig\tchrMapping\n");
    for (BasicFeature gene : genes) {
      geneContigWriter.write(
          gene.id + "\t" + gene.contig + "\t" + seqIdTochrMapping.get(gene.contig) + "\n");
    }
    geneContigWriter.close();

    FileWriter geneIdMappingWriter = new FileWriter(geneIdMappingFile);
    geneIdMappingWriter.write("id\txRefGeneId\tonMainContig\n");
    for (BasicFeature gene : genes) {
      geneIdMappingWriter.write(gene.id + "\t" + gene.xRefGeneId + "\t" + gene.onMainContig + "\n");
    }
    geneIdMappingWriter.close();
  }
}
