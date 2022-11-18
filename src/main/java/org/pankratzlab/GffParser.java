package org.pankratzlab;

import htsjdk.samtools.util.LineReader;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.tribble.readers.LineIterator;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

public class GffParser {
  private final ArrayList<Gff3Feature> features= new ArrayList<Gff3Feature>(1_000_000);
  public GffParser(String filename) {
    File inputFile = new File(filename);
    if (!inputFile.exists()) {
      throw new IllegalArgumentException("Input file does not exist");
    }
    Gff3Codec gff3Codec = new Gff3Codec();
    if (!gff3Codec.canDecode(filename)) {
      throw new IllegalArgumentException("Gff3Codec says it cannot decode this file!");
    }

    try (InputStream instream = new FileInputStream(inputFile)) {
      LineIterator iterator = new AsciiLineReaderIterator(AsciiLineReader.from(instream));
      final AbstractFeatureReader<Gff3Feature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(filename, gff3Codec, false);
      for(Gff3Feature f : reader.iterator()) {
        features.add(f);
      }
      System.out.println("something?");
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
