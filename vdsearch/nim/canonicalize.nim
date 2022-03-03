import bioseq
import cello
import strutils

import nimpy


proc minimalCanonicalRotation*(x:  Dna): Dna = 
  var rotatedCanonical: Dna
  result = x # we'll start by assuming the assuming the sequence is already rotationally canonical
  var rc = x.reverseComplement
  for i in 0..x.high:
    rotatedCanonical = min($rotate($x, i), $rotate($rc, i)).toDna
    if rotatedCanonical < result:
      result = rotatedCanonical

proc canonicalize*(infile: string, outfile: string) {.exportpy.}=
  let outfileFile = open(outfile, fmWrite) # the output file as an opend File object
  defer: outfileFile.close()
  for record in readFasta[Dna](infile):
    writeLine outfileFile, toRecord[Dna](record.sequence.minimalCanonicalRotation, record.description).asFasta