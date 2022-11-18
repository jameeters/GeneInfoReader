package org.pankratzlab;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Constants;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class GffParser {
  private final ArrayList<Gff3Feature> genes = new ArrayList<>(1_000_000);
  public GffParser(String filename) {
    File inputFile = new File(filename);
    if (!inputFile.exists()) {
      throw new IllegalArgumentException("Input file does not exist");
    }

    // htsjdk will not let us throw these away
    Set<String> requiredAttributes = new HashSet<>(Arrays.asList(
        Gff3Constants.ID_ATTRIBUTE_KEY,
        Gff3Constants.PARENT_ATTRIBUTE_KEY,
        Gff3Constants.NAME_ATTRIBUTE_KEY
    ));

    Gff3Codec gff3Codec = new Gff3Codec(Gff3Codec.DecodeDepth.DEEP, S -> !requiredAttributes.contains(S));
    if (!gff3Codec.canDecode(filename)) {
      throw new IllegalArgumentException("Gff3Codec says it cannot decode this file!");
    }

    Set<String> typesWeCareAbout = new HashSet<>(List.of("gene"));

    try {
      final AbstractFeatureReader<Gff3Feature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(filename, gff3Codec, false);
      reader.iterator().stream()
          .filter(f -> typesWeCareAbout.contains(f.getType()))
          .forEach(genes::add);

      System.out.println("Finished loading " + genes.size() + " genes");
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
