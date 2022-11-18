package org.pankratzlab;

public class Main {
  public static void main(String[] args) {
    String gffFile = args[0];
    GffParser parser = new GffParser(gffFile);
  }
}