package org.pankratzlab;

import java.io.IOException;
import java.nio.file.Path;

// https://ftp.ncbi.nlm.nih.gov//genomes/all/annotation_releases/9606/109.20210514/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
public class Main {
  public static void main(String[] args) {
    Path gffFile = Path.of(args[0]);
    Path outputDir = Path.of("/tmp");
    boolean qc = true;
    boolean bedExons = false;
    boolean bedIntrons = false;
    boolean bedAll = false;
    if (args.length >= 2) {
      for(String a : args) {
        if (a.equals("-noqc")) {
          qc = false;
        } else if (a.equals("-bedExons")) {
          bedExons = true;
        } else if(a.equals("-bedIntrons")) {
          bedIntrons = true;
        } else if (a.startsWith("-out=")) {
          outputDir = Path.of(a.replace("-out=", ""));
        } else if (a.equals("-bedAll")) {
          bedAll = true;
        }
      }
    }

    Aggregator aggregator = new Aggregator(gffFile, outputDir);

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
    } else {
      aggregator.writeSerializedGeneTrack();
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