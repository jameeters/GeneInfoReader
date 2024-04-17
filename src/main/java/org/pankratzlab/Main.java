package org.pankratzlab;

import java.io.IOException;

public class Main {
  public static void main(String[] args) {
    String gffFile = args[0];
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
        }
      }
    }

    Aggregator aggregator = new Aggregator(gffFile);

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