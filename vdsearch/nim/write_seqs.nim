import bioseq
import nimpy
import std/[sets, tables, os, strutils]

proc write_seqs*(infile: string, outfile: string, ids: seq[string]): void {.exportpy.} =
  let idSet = toHashSet(ids)
  let outfileFile = open(outfile, fmWrite) # the output file as an opend File object
  defer: outfileFile.close()
  for record in readFasta[Dna](infile):
    # Infernal only reports the ID, not the full header so we have to parse it
    # I'm not a fan of this trying to guess FASTA header format
    if record.description.splitWhitespace()[0] in idSet:
      outfileFile.writeLine(record.asFasta())

proc write_clusters*(infile: string, outdir: string, mapping: Table[string, string], clusters: seq[string]): void {.exportpy.} =
  ## Write out a fasta file for each cluster.
  ## 
  ## Note: `outdir` must exist and be writable or the function will fail.  

  var outfiles = initTable[string, File]()

  try:
    for cluster in clusters:
      outfiles[cluster] = open(outdir / cluster & ".fasta", fmWrite) # the output file as an opend File object

    for record in readFasta[Dna](infile):
      let clusterFile = outfiles[mapping[record.description]]
      clusterFile.writeLine(record.asFasta())
  
  finally:
    for cluster in clusters:
      outfiles[cluster].close()
