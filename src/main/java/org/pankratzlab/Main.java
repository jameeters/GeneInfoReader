package org.pankratzlab;

import java.io.IOException;
import java.nio.file.Path;

public class Main {
  public static void main(String[] args) {
    Path gffFile = Path.of(args[0]);
    Path outputDir = Path.of("/tmp");
    boolean qc = true;
    boolean bedExons = false;
    boolean bedIntrons = false;
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
        }
      }
    }

    Aggregator aggregator = new Aggregator(gffFile, outputDir);

    aggregator.findGenesAndExons();
    if (bedIntrons) {
      aggregator.findGenesAndIntrons();
    }
    aggregator.computeXRefMap();

    if (bedExons || bedIntrons) {
      aggregator.writeBedFile();
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