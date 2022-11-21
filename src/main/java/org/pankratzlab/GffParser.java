package org.pankratzlab;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class GffParser {
  Aggregator aggregator = new Aggregator();

  public GffParser(String filename) {
    File inputFile = new File(filename);
    if (!inputFile.exists()) {
      throw new IllegalArgumentException("Input file does not exist");
    }

    Gff3Codec gff3Codec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);
    if (!gff3Codec.canDecode(filename)) {
      throw new IllegalArgumentException("Gff3Codec says it cannot decode this file!");
    }

    try {
      final AbstractFeatureReader<Gff3Feature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(filename, gff3Codec, false);
      reader.iterator().stream()
          .map(Gff3Feature::getBaseData)
          .filter(f -> !f.getAttribute("pseudo").equals(List.of("true")))
          .forEach(aggregator::add);

      System.out.println("Finished loading " + aggregator.featureMap.size() + " features");
      reader.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
