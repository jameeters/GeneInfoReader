package org.pankratzlab;

import java.io.IOException;

public class Main {
  public static void main(String[] args) {
    String gffFile = args[0];
    boolean qc = true;
    if(args.length == 2){
      if(args[1].equals("-noqc")){
        qc = false;
      }
    }
    GffParser parser = new GffParser(gffFile);
    parser.aggregator.findGenesAndExons();
    parser.aggregator.computeXRefMap();

    if(qc) {
      try {
        parser.aggregator.writeTestOutput();
        parser.aggregator.writeQcOutput();
      } catch (IOException e) {
        e.printStackTrace();
      }
    }

    System.out.println("done");
  }
}