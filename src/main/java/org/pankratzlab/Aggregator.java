package org.pankratzlab;

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

import org.pankratzlab.common.filesys.GeneData;
import org.pankratzlab.common.filesys.GeneSet;
import org.pankratzlab.common.filesys.GeneTrack;

import htsjdk.tribble.gff.Gff3BaseData;
import htsjdk.tribble.gff.Gff3Feature;

public class Aggregator {
  final GffParser parser;

  final Map<String, BasicFeature> featureMap = new HashMap<>();
  final Set<BasicFeature> genes = new HashSet<>();
  final Set<String> duplicateIds = new HashSet<>();

  final static String BAD_OR_MISSING = "BAD_OR_MISSING";
  Map<String, GeneGrouping> geneGroupingsByXRefGeneId = new TreeMap<>();

  public Aggregator(String gffFilename) {
    this.parser = new GffParser(gffFilename, this::add);
    System.out.println("Finished loading " + featureMap.size() + " features");
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
    System.out.println("Finding genes and exons...");
    for (BasicFeature feat : featureMap.values()) {
      if (feat.isGene()) {
        this.genes.add(feat);
        feat.getDescendantExons();
      }
    }
  }

  public void computeXRefMap() {
    System.out.println("Computing gene groups based on xRefGeneId...");
    for (BasicFeature gene : genes) {
      geneGroupingsByXRefGeneId.computeIfAbsent(gene.xRefGeneId,
                                                (String xRefId) -> new GeneGrouping(xRefId,
                                                                                    genes.parallelStream()
                                                                                         .filter(g -> g.xRefGeneId.equals(xRefId))
                                                                                         .collect(Collectors.toSet()))

      );
      if (geneGroupingsByXRefGeneId.size() % 2000 == 0) {
        System.out.println(geneGroupingsByXRefGeneId.size() + " groups computed");
      }
    }
  }

  public void writeSerializedGeneTrack() {
    String geneSetFile = "/tmp/geneset.ser";
    String geneTrackFile = "/tmp/GeneTrack.ser";

    System.out.println("Creating GeneTrack...");
    List<GeneData> geneDatas = geneGroupingsByXRefGeneId.values().stream()
                                                        .map(GeneGrouping::getMainContigGenes)
                                                        .flatMap((Set<BasicFeature> s) -> s.stream()
                                                                                           .map(BasicFeature::toGeneData))
                                                        .collect(Collectors.toList());

    GeneSet geneSet = new GeneSet(geneDatas);
    geneSet.serialize(geneSetFile);

    GeneTrack geneTrack = new GeneTrack(geneSetFile);
    geneTrack.serialize(geneTrackFile);
  }

  void writeQcOutput() throws IOException {
    System.out.println("Writing QC files...");
    String genesAndExonsFile = "/tmp/geneinfo";
    String chrGeneCountsFile = "/tmp/chrGeneCounts.tsv";
    String seqIdCountsFile = "/tmp/seqIdCounts.tsv";
    String duplicateIdsFile = "/tmp/duplicateIds.tsv";
    String genesContigsFile = "/tmp/genesContigs.tsv";
    String geneIdMappingFile = "/tmp/geneIdMapping.tsv";
    String geneGroupingsFile = "/tmp/geneGroupings.tsv";

    List<String> fileNames = List.of(genesAndExonsFile, chrGeneCountsFile, seqIdCountsFile,
                                     duplicateIdsFile, genesContigsFile, geneIdMappingFile,
                                     geneGroupingsFile);

    for (String fileName : fileNames) {
      File f = new File(fileName);
      boolean result = Files.deleteIfExists(f.toPath());
      if (result) {
        System.out.println("deleted " + fileName);
      }
    }

    FileWriter genesAndExonsWriter = new FileWriter(genesAndExonsFile);
    for (GeneGrouping gg : geneGroupingsByXRefGeneId.values()) {
      for (BasicFeature gene : gg.getMainContigGenes()) {
        genesAndExonsWriter.write(gene.name + "  " + gene.start + "  " + gene.end + "\n");
        for (BasicFeature exon : gene.getDescendantExons()) {
          genesAndExonsWriter.write("   |" + exon.id + "  " + exon.start + "  " + exon.end + "\n");
        }
      }
    }
    genesAndExonsWriter.close();

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
      geneContigWriter.write(gene.id + "\t" + gene.contig + "\t"
                             + seqIdTochrMapping.get(gene.contig) + "\n");
    }
    geneContigWriter.close();

    FileWriter geneIdMappingWriter = new FileWriter(geneIdMappingFile);
    geneIdMappingWriter.write("id\txRefGeneId\tonMainContig\n");
    for (BasicFeature gene : genes) {
      geneIdMappingWriter.write(gene.id + "\t" + gene.xRefGeneId + "\t" + gene.onMainContig + "\n");
    }
    geneIdMappingWriter.close();

    FileWriter geneGroupingsWriter = new FileWriter(geneGroupingsFile);
    geneGroupingsWriter.write("xRefGeneId\ttotalGenes\tmainContigGenes\n");
    for (GeneGrouping gg : this.geneGroupingsByXRefGeneId.values()) {
      geneGroupingsWriter.write(gg.geneId + "\t" + gg.countTotalGenes() + "\t"
                                + gg.countMainContigGenes() + "\n");
    }
    geneGroupingsWriter.close();
  }
}
