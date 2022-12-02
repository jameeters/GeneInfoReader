package org.pankratzlab;

import java.io.IOException;

public class Main {
  public static void main(String[] args) {
    String gffFile = args[0];
    boolean qc = true;
    if (args.length == 2) {
      if (args[1].equals("-noqc")) {
        qc = false;
      }
    }

    Aggregator aggregator = new Aggregator(gffFile);

    aggregator.findGenesAndExons();
    aggregator.computeXRefMap();

    aggregator.writeSerializedGeneTrack();

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