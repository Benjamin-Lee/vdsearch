import nimpy

import bioseq
import strutils


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
  var mers = 2
  while monomerized != Dna"":
      newMonomer = monomerized.cirit(seedLen, minIdentity)
      if newMonomer == Dna"":
        break
      monomerized = newMonomer
      inc mers

  return toRecord[Dna](monomerized, record.description)

proc find_circs*(infile: string, outfile: string, seedLen: int = 10, minIdentity: float = 0.95, reportMultimers: bool = false): int {.exportpy.} =
    let outfileFile = open(outfile, fmWrite) # the output file as an opend File object
    var monomerized: Record[Dna]
    var count = 0
    for record in readFasta[Dna](infile):
      monomerized = record.monomerize(seedLen, minIdentity)
      if monomerized.len > 0:
        writeLine(outfileFile, monomerized.asFasta())
        inc count
    return count