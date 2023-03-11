import strutils
import strformat
import hashes
import tables
import os

type 
  Monomer* = distinct char
  Prot* = distinct string
  Dna* = distinct string
  Rna* = distinct string
  BioString* = Dna|Rna|Prot
  NucleicAcid = Dna|Rna
  Record*[T: BioString] = ref object
    sequence*: T
    description*: string
    quality*: string
  
proc `==`*(a, b: Monomer): bool {.borrow.}
proc `$`*(a: Monomer): string {.borrow.}

template defineStrOprs(typ: typedesc) {.dirty.} =
  proc `$`*(x: typ): string {.borrow.}
  template `&`*(x, y: typ): typ = typ(x.string & y.string)
  proc `&=`*(x: var typ, y: typ) {.borrow.}
  template `[]`*[T, U](x: typ; h: HSlice[T, U]): typ = typ(x.string[h])
  template `[]`*(x: typ, i: int): Monomer = Monomer(x.string[i])
  template `[]`*(x: typ; i: BackwardsIndex): Monomer = Monomer(x.string[i])
  template `[]=`*[T, U](s: var typ, x: HSlice[T, U], b: typ) = s.string[x] = b.string
  template `[]=`*(s: typ; i: int; val: Monomer) = s.string[i] = val.char
  proc `==`*(x, y: typ): bool {.borrow.}
  proc `<`*(x, y: typ): bool {.borrow.}
  proc `<=`*(x, y: typ): bool {.borrow.}
  proc high*(x: typ): int {.borrow.}
  proc low*(x: typ): int {.borrow.}
  proc len*(x: typ): int {.borrow.}
  proc hash*(x: typ): Hash {.borrow.}
  proc startsWith*(x, y: typ): bool {.borrow.}
  proc endsWith*(x, y: typ): bool {.borrow.}
  proc continuesWith*(x, y: typ, start: Natural): bool {.borrow.}
  proc contains*(x, y: typ): bool {.borrow.}
  converter toSequence*(x: Record[typ]): typ = x.sequence
  iterator items*(x: typ): Monomer {.inline.}  =
    for base in x.string:
      yield Monomer(base)
  iterator pairs*(x: typ): tuple[key: int, val: Monomer] {.inline.} =
    for i, base in x.string:
      yield (i, x[i])
  
defineStrOprs(Dna)
defineStrOprs(Rna)
defineStrOprs(Prot)

proc toDna*(x: string): Dna = x.Dna
proc toRna*(x: string): Rna = x.Rna
proc toProtein*(x: string): Prot = x.Prot
proc toRecord*[T: BioString](x: T, description = "", quality = ""): Record[T] = 
  Record[T](sequence: x, description: description, quality: quality)
proc newRecord*[T: BioString](x: string, description="", quality=""): Record[T] = 
  Record[T](sequence: x.T, description: description, quality: quality)
  
proc hash*(x: Record): Hash =
  var h: Hash = 0
  h = h !& hash(x.sequence)
  h = h !& hash(x.description)
  h = h !& hash(x.quality)
  result = !$h
proc `==`*(x, y: Record): bool = 
  x.sequence == y.sequence and x.description == y.description and x.quality == y.quality  

proc asFasta*(x: Record): string = 
  result &= ">" & x.description & "\n"
  result &= x.sequence.string
proc asFastq*[T: Dna|Rna](x: Record[T]): string = 
  if x.sequence.len != x.quality.len:
    raise newException(CatchableError, "Quality string must be the same length as sequence")
  result &= "@" & x.description & "\n"
  result &= x.sequence.string & "\n"
  result &= "+\n"
  result &= x.quality
    
proc gcContent*(x: NucleicAcid): float =
  for letter in $x:
    case letter:
      of 'A', 'T', 'U', 'a', 't', 'u':
        continue
      of 'G', 'C', 'g', 'c':
        result += 1
      else: 
        continue
  result /= x.len.float

iterator kmers*[T: BioString](x: T, k: Positive): T =
  for i in 0..(x.len - k):
    yield x[i ..< i + k]
template canonical*[T: Dna|Rna](x: T): T = 
  min(x, x.reverseComplement)
iterator canonicalKmers*[T: BioString](x: T, k: Positive): T =
  for kmer in x.kmers(k):
    yield canonical(kmer)
proc countKmers*[T: Dna|Rna](x: T, k: Positive): CountTable[T] =
  for kmer in x.kmers(k):
    result.inc(kmer)
func totalKmers*[T: Dna|Rna](x: T, k: Positive): Natural =
  if k > x.len:
    return 0
  else:
    return x.len - k + 1

proc reverseComplement*[T: Dna|Rna](x: T): T =
  result = newString(x.len).T
  var i = x.high
  var j = 0
  while i >= 0:
    result.string[i] = case x.string[j]:
      of 'A':
        when T is Dna: 'T' else: 'U'
      of 'T', 'U':
        'A'
      of 'G':
        'C'
      of 'C':
        'G'
      of 'Y':
        'R'
      of 'R':
        'Y'
      of 'S':
        'S'
      of 'W':
        'W'
      of 'K':
        'M'
      of 'M':
        'K'
      of 'B':
        'V'
      of 'D':
        'H'
      of 'H':
        'D'
      of 'V':
        'B'
      of 'N':
        'N'
      else:
        raise newException(CatchableError, "Invalid character in sequence")
    inc(j)
    dec(i)
template `~`*[T: Dna|Rna](x: T): T = x.reverseComplement

proc overlapsWithStartOf*(a: Dna, b: Dna, k: int): bool {.inline.}= 
  ## Checks whether the end of `a` overlaps with the start (left) of `b` by at least `k` bases.
  runnableExamples:
    import bioseq
    # the last three of `a` and first three of `b` overlap
    assert "ATGCAGA".toDna.overlapsWithStartOf("AGATTAGATA".toDna, 3) == true
    # here, the last three are not equivalent to the first three but there is still an overlap
    assert "GGCCAAGCCC".toDna.overlapsWithStartOf("GCCCAGGTATGC".toDna, 3) == true
    # note that CGA (first three bases of `b`) also occurs at the start of `a`
    assert "CGATATTTCGATA".toDna.overlapsWithStartOf("GATATCAGAA".toDna, 3) == true
    # if `b` is a substring of `a` it doesn't count as overlapping
    assert "CCATG".toDna.overlapsWithStartOf("CATG".toDna, 3) == false
  var posInB: int # where we are in `b`
  var matchFound: bool # whether we're exiting because we found a match
  for i in countdown(a.len - k, a.low):
    # return early if you've reached the end of possible overlaps
    if a.len - i >= b.len:
      return false
    # reset the iteration variables
    posInB = 0
    matchFound = true
    # starting from i, iterate to the end of the sequence
    for posInA in i..a.high:
      # if the bases differ, it's no match
      if a[posInA] != b[posInB]:
        matchFound = false
        break
      # otherwise, check the next base
      inc(posInB)
    # exit early
    if matchFound:
      return true
  return false

iterator readFasta*[T: BioString](filename: string|File): Record[T] =
  ## Iterate over the lines in a FASTA file, yielding one record at a time 
  var description = ""
  var sequence = ""
  for line in filename.lines:
    if line.startsWith(">"):
      if sequence != "":
        yield newRecord[T](sequence, description=description)
        sequence = ""
      description = line[1..line.high]
    else:
      sequence &= line.strip
  yield newRecord[T](sequence, description=description)
  sequence = ""
  
iterator readFastq*[T: NucleicAcid](filename: string|File): Record[T] =
  var linesRead = 0
  var description = ""
  var sequence = ""
  for line in filename.lines:
    if linesRead mod 4 == 0:
      description = line
      description.delete(0..0) # remove the @
    elif linesRead mod 4 == 1:
      sequence = line
    elif linesRead mod 4 == 3:
      yield newRecord[T](sequence, description=description, quality=line)
      description = ""
      sequence = ""
    linesRead += 1

iterator readFastx*[T: NucleicAcid](path: string): Record[T] =
  # Use the correct parser depending on the file extension
  var mode = if path.splitFile.ext.toLowerAscii in [".fastq", ".fq"]: "fastq" else: "fasta"
  if mode == "fasta":
    for record in readFasta[T](path): yield record
  elif mode == "fastq":
    for record in readFastq[T](path): yield record
