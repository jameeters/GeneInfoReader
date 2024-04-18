package org.pankratzlab;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.StringJoiner;
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
  private boolean genesFound = false;
  final Set<String> duplicateIds = new HashSet<>();
  final Path outputDir;
  final static String BAD_OR_MISSING = "BAD_OR_MISSING";
  Map<String, GeneGrouping> geneGroupingsByXRefGeneId = new TreeMap<>();

  public Aggregator(Path gffFilename, Path outputDir) {
    this.outputDir = outputDir;
    this.parser = new GffParser(gffFilename.toString(), this::add);
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

  private void findGenes() {
    System.out.println("Finding genes...");
    for (BasicFeature feat : featureMap.values()) {
      if (feat.isGene()) {
        this.genes.add(feat);
      }
    }
    this.genesFound = true;
  }

  void findGenesAndExons() {
    if (!genesFound) {
      this.findGenes();
    }
    System.out.println("Finding exons...");
    this.genes.forEach(BasicFeature::getDescendantExons);
  }

  void findGenesAndIntrons() {
    if (!genesFound) {
      this.findGenes();
    }
    System.out.println("Finding introns...");
    this.genes.forEach(BasicFeature::getDescendantIntrons);
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
    Path geneSetFile = outputDir.resolve("geneset.ser");
    Path geneTrackFile = outputDir.resolve("GeneTrack.ser");

    System.out.println("Creating GeneTrack...");
    List<GeneData> geneDatas = geneGroupingsByXRefGeneId.values().stream()
                                                        .map(GeneGrouping::getMainContigGenes)
                                                        .flatMap((Set<BasicFeature> s) -> s.stream()
                                                                                           .map(BasicFeature::toGeneData))
                                                        .collect(Collectors.toList());

    GeneSet geneSet = new GeneSet(geneDatas);
    geneSet.serialize(geneSetFile.toString());

    GeneTrack geneTrack = new GeneTrack(geneSetFile.toString());
    geneTrack.serialize(geneTrackFile.toString());
  }

  public void writeBedFile(boolean includeExons, boolean includeIntrons) {
    StringJoiner filename = new StringJoiner("_");
    if (includeExons) filename.add("exons");
    if (includeIntrons) filename.add("introns");
    Path bedFile = outputDir.resolve(filename + ".bed");
    PrintWriter writer = org.pankratzlab.common.Files.getAppropriateWriter(bedFile.toString());
    geneGroupingsByXRefGeneId.values().stream().sorted(GeneGrouping::compareTo)
        .forEachOrdered((GeneGrouping geneGrouping) -> {
          try {
            BasicFeature mainContigGene = geneGrouping.getMainContigGenes().stream().findFirst()
                .get();
            mainContigGene.toBedLines(includeExons, includeIntrons).forEach(writer::println);
          } catch (NoSuchElementException e) {
            System.out.println("No main contig gene found for group " + geneGrouping.geneId);
          }
        });
    writer.close();
  }

  void writeQcOutput() throws IOException {
    System.out.println("Writing QC files...");
    Path genesAndExonsFile = outputDir.resolve("geneinfo");
    Path chrGeneCountsFile = outputDir.resolve("chrGeneCounts.tsv");
    Path seqIdCountsFile = outputDir.resolve("seqIdCounts.tsv");
    Path duplicateIdsFile = outputDir.resolve("duplicateIds.tsv");
    Path genesContigsFile = outputDir.resolve("genesContigs.tsv");
    Path geneIdMappingFile = outputDir.resolve("geneIdMapping.tsv");
    Path geneGroupingsFile = outputDir.resolve("geneGroupings.tsv");

    List<Path> fileNames = List.of(genesAndExonsFile, chrGeneCountsFile, seqIdCountsFile,
                                     duplicateIdsFile, genesContigsFile, geneIdMappingFile,
                                     geneGroupingsFile);

    for (Path fileName : fileNames) {
      File f = fileName.toFile();
      boolean result = Files.deleteIfExists(f.toPath());
      if (result) {
        System.out.println("deleted " + fileName);
      }
    }

    FileWriter genesAndExonsWriter = new FileWriter(genesAndExonsFile.toFile());
    for (GeneGrouping gg : geneGroupingsByXRefGeneId.values()) {
      for (BasicFeature gene : gg.getMainContigGenes()) {
        genesAndExonsWriter.write(gene.name + "  " + gene.start + "  " + gene.end + "\n");
        for (BasicFeature exon : gene.getDescendantExons()) {
          genesAndExonsWriter.write("   |" + exon.id + "  " + exon.start + "  " + exon.end + "\n");
        }
      }
    }
    genesAndExonsWriter.close();

    FileWriter geneCountsWriter = new FileWriter(chrGeneCountsFile.toFile());
    FileWriter seqIdCountsWriter = new FileWriter(seqIdCountsFile.toFile());

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

    FileWriter duplicateIdsWriter = new FileWriter(duplicateIdsFile.toFile());
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

    FileWriter geneContigWriter = new FileWriter(genesContigsFile.toFile());
    geneContigWriter.write("id\tcontig\tchrMapping\n");
    for (BasicFeature gene : genes) {
      geneContigWriter.write(gene.id + "\t" + gene.contig + "\t"
                             + seqIdTochrMapping.get(gene.contig) + "\n");
    }
    geneContigWriter.close();

    FileWriter geneIdMappingWriter = new FileWriter(geneIdMappingFile.toFile());
    geneIdMappingWriter.write("id\txRefGeneId\tonMainContig\n");
    for (BasicFeature gene : genes) {
      geneIdMappingWriter.write(gene.id + "\t" + gene.xRefGeneId + "\t" + gene.onMainContig + "\n");
    }
    geneIdMappingWriter.close();

    FileWriter geneGroupingsWriter = new FileWriter(geneGroupingsFile.toFile());
    geneGroupingsWriter.write("xRefGeneId\ttotalGenes\tmainContigGenes\n");
    for (GeneGrouping gg : this.geneGroupingsByXRefGeneId.values()) {
      geneGroupingsWriter.write(gg.geneId + "\t" + gg.countTotalGenes() + "\t"
                                + gg.countMainContigGenes() + "\n");
    }
    geneGroupingsWriter.close();
  }
}
