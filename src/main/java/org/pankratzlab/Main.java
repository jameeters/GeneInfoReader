package org.pankratzlab;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

public class Main {
  public static void main(String[] args) {
    // @formatter:off
    String usage = "\n" + "GeneInfoReaderUsage: \n"
                   + "dir=./  specify working directory (optional)\n"
                   + "inputFile=... specify gff3 input file (required)\n"
                   + "qc=true whether to output qc files (optional)\n" + "\n";
    // @formatter:on

    Path inputFile = null;
    boolean qc = true;
    Path dir = Paths.get("./");

    for (String arg : args) {
      if (arg.startsWith("dir=")) {
        dir = Paths.get(arg.split("=")[1]);
      } else if (arg.startsWith("inputFile=")) {
        Path p = Paths.get(arg.split("=")[1]);
        inputFile = dir.resolve(p);
      } else if (arg.startsWith("qc=")) {
        qc = Boolean.parseBoolean(arg.split("=")[1]);
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

    Aggregator aggregator = new Aggregator(dir, inputFile);

    aggregator.findGenesAndExons();
    aggregator.computeXRefMap();

    aggregator.writeSerializedGeneTrack();

    try {
      aggregator.writeGenesXlnFile();
    } catch (IOException e) {
      e.printStackTrace();
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