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

proc canonicalize*(infile: string, outfile: string, minLen: Natural = 1, maxLen: Natural = high(int)) {.exportpy.}=

  # Warn the user if they compiled wrong
  if not defined(danger):
    echo "Not compiled with -d:danger. This will likely cause severe slowdowns."

  let outfileFile = open(outfile, fmWrite) # the output file as an opend File object
  defer: outfileFile.close()
  var tmp = toRecord[Dna](Dna"", "") # only allocate a temporary record once
  for record in readFasta[Dna](infile):
    
    # bail to prevent slowdowns on giant sequences
    if record.len < minLen or record.len > maxLen:
      continue

    tmp.sequence = record.sequence.minimalCanonicalRotation
    tmp.description = record.description
    writeLine outfileFile, tmp.asFasta