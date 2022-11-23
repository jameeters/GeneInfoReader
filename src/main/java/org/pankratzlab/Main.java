package org.pankratzlab;

import java.io.IOException;

public class Main {
  public static void main(String[] args) {
    String gffFile = args[0];
    GffParser parser = new GffParser(gffFile);
    parser.aggregator.findGenesAndExons();
    try {
      parser.aggregator.writeTestOutput();
      parser.aggregator.writeQcOutput();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}