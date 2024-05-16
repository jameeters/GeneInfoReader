package org.pankratzlab;

import java.io.IOException;
import java.nio.file.Path;

// https://ftp.ncbi.nlm.nih.gov//genomes/all/annotation_releases/9606/109.20210514/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
public class Main {
  public static void main(String[] args) {
    // @formatter:off
    String usage = "\n" + "GeneInfoReaderUsage: \n"
                   + "inputFile=... specify gff3 input file (required)\n"
                   + "out=/tmp/ specify an output directory (optional)\n"
                   + "-noqc provide this flag to suppress QC output\n"
                   + "-geneTrack provide this flag to generate a serialized GeneTrack for "
                   + "Genvisis\n"
                   + "-genesXln provide this flag to generate an xln file of genes\n"
                   + "-bedExons provide this flag to generate a bed file of exons\n"
                   + "-bedIntrons provide this flag to generate a bed file of introns\n"
                   + "-bedAll provide this flag to generate three bed files. One of exons, one of "
                   + "introns, and one containing both. \n" + "\n";
    // @formatter:on

    Path inputFile = null;
    Path outputDir = Path.of("/tmp");
    boolean qc = true;
    boolean geneTrack = false;
    boolean genesXln = false;
    boolean bedExons = false;
    boolean bedIntrons = false;
    boolean bedAll = false;

    for (String arg : args) {
      if (arg.startsWith("inputFile=")) {
        inputFile = Path.of(arg.split("=")[1]);
      } else if (arg.equals("qc=")) {
        qc = Boolean.parseBoolean(arg.split("=")[1]);
      } else if (arg.equals("-noqc")) {
        qc = false;
      } else if (arg.equals("-geneTrack")) {
        geneTrack = true;
      } else if (arg.equals("-genesXln")) {
        genesXln = true;
      } else if (arg.equals("-bedExons")) {
        bedExons = true;
      } else if (arg.equals("-bedIntrons")) {
        bedIntrons = true;
      } else if (arg.startsWith("out=")) {
        outputDir = Path.of(arg.replace("out=", ""));
      } else if (arg.equals("-bedAll")) {
        bedAll = true;
      } else {
        System.out.println(usage);
        System.exit(1);
      }
    }
    if (inputFile == null) {
      System.err.println("No input file provided!");
      System.out.println(usage);
      System.exit(1);
    }

    if (!(geneTrack || genesXln || bedIntrons || bedExons || bedAll)) {
      System.out.println("You haven't asked for any output...");
      System.out.println(usage);
      System.exit(0);
    }

    Aggregator aggregator = new Aggregator(inputFile, outputDir);

    aggregator.findGenesAndExons();
    if (bedIntrons || bedAll) {
      aggregator.findGenesAndIntrons();
    }
    aggregator.computeXRefMap();

    if (bedAll) {
      aggregator.writeBedFile(false, true);
      aggregator.writeBedFile(true, false);
      aggregator.writeBedFile(true, true);
    } else if (bedExons || bedIntrons) {
      aggregator.writeBedFile(bedExons, bedIntrons);
    }

    if (geneTrack) {
      aggregator.writeSerializedGeneTrack();
    }

    if (genesXln) {
      try {
        aggregator.writeGenesXlnFile();
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
    if (qc) {
      try {
        aggregator.writeQcOutput();
      } catch (IOException e) {
        e.printStackTrace();
      }
    }

    System.out.println("done");
  }
}