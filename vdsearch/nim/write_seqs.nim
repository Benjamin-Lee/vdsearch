import bioseq
import nimpy
import std/sets

proc write_seqs*(infile: string, outfile: string, ids: seq[string]): void {.exportpy.} =
  let idSet = toHashSet(ids)
  let outfileFile = open(outfile, fmWrite) # the output file as an opend File object
  for record in readFasta[Dna](infile):
    if record.description in idSet:
      writeLine(outfileFile, record.asFasta())