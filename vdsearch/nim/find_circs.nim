import nimpy
import std/[os, monotimes, times]
import bioseq
import strutils
from canonicalize import minimalCanonicalRotation


func cirit(x: Dna, seedLen = 10, minIdentity=0.95): Dna = 

  # prevent an error from the seed being larger than the sequence
  if x.len < seedLen:
    return Dna""
  
  # find the seed (short sequence at the end that must match exactly)
  let seed = x[^seedLen..^1]

  # search the sequence for it
  let idx = ($x).find($seed, last = x.len - seedLen)
  
  # find returns -1 if not found, so if it is -1, there's no match
  if idx == -1:
    return Dna""
  
  # compare the start of the sequence to the seed
  let maxMismatches = (idx.float * (1.0 - minIdentity)).int
  var mismatches = 0
  for i in 0..idx:
    if x[i] != x[ x.high - idx + i - seedLen + 1]:
      inc mismatches
    if mismatches > maxMismatches:
      return Dna""

  return x[0..x.high - idx - seedLen]

func monomerize(record: Record[Dna], seedLen=10, minIdentity=0.95): Record[Dna] =
  # first, try to see if the sequence is a circRNA
  var monomerized = record.cirit(seedLen, minIdentity)

  # if not, go on
  if monomerized == Dna"":
    return toRecord[Dna](monomerized, record.description)

  var newMonomer: Dna
  while monomerized != Dna"":
    newMonomer = monomerized.cirit(seedLen, minIdentity)
    if newMonomer == Dna"":
      break
    monomerized = newMonomer

  return toRecord[Dna](monomerized, record.description)

proc find_circs*(infile: string,
                 outfile: string,
                 seedLen: Natural = 10,
                 minIdentity: float = 0.95,
                 canonicalize: bool = true,
                 outTsv: bool = true,
                 minLen: Natural = 1,
                 maxLen: Natural = high(int),
                ): (int, int, int) {.exportpy.} =

  let outfileFile = open(outfile, fmWrite) # the output file as an opend File object
  defer: outfileFile.close()

  var outTsvFile: File
  defer: outTsvFile.close()

  if outTsv:
    outTsvFile = open(changeFileExt(outfile, "tsv"), fmWrite) # the output file as an opend File object
    outTsvFile.writeLine("id", "\t", "ratio", "\t", "original_len", "\t", "monomerized_len")

  var monomerized: Record[Dna]
  var count = 0
  var total = 0
  var originalLen = 0
  let startTime = getMonoTime()

  for record in readFasta[Dna](infile):
    originalLen = record.len
    monomerized = record.monomerize(seedLen, minIdentity)

    if canonicalize:
      monomerized = toRecord[Dna](minimalCanonicalRotation(monomerized), monomerized.description)

    if minLen <= monomerized.len and monomerized.len <= maxLen:
      writeLine(outfileFile, monomerized.asFasta())
      inc count
      
      # Write out the ratio between the original and monomerized sequence length to a TSV file
      # This is useful since finding the original might take a long time
      if outTsv:
        writeLine(outTsvFile, monomerized.description, "\t", originalLen / monomerized.len, "\t", originalLen, "\t", monomerized.len)
    
    inc total

  return (count: count, total: total, duration: (getMonoTime() - startTime).inMilliseconds.int)