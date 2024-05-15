package org.pankratzlab;

import java.io.File;
import java.io.IOException;
import java.util.function.Consumer;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;

public class GffParser {

  public GffParser(String filename, Consumer<Gff3Feature> featureConsumer) {
    File inputFile = new File(filename);
    if (!inputFile.exists()) {
      throw new IllegalArgumentException("Input file does not exist: " + filename);
    }

    Gff3Codec gff3Codec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);
    if (!gff3Codec.canDecode(filename)) {
      throw new IllegalArgumentException("Gff3Codec says it cannot decode this file!");
    }

    try {
      final AbstractFeatureReader<Gff3Feature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(filename,
                                                                                                             gff3Codec,
                                                                                                             false);
      reader.iterator().stream().forEach(featureConsumer);

      reader.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
