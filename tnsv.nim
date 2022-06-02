import argparse
import hts/vcf
import hts/private/hts_concat
import strformat
import sets
import std/random
import strutils
import lapper
import tables

randomize(42)

type expanded_sv = object
  start: int
  stop: int

proc start(e:expanded_sv): int {.inline.} = return e.start
proc stop(e:expanded_sv): int {.inline.} = return e.stop

proc key(v:Variant): string =
  return &"{v.CHROM}/{v.start}/{v.REF}/{v.ALT[0]}"

template stripChr*[T:string|cstring](s:T): string =
  if s.len > 3 and ($s).startswith("chr"): ($s)[3..<s.len] else: $s

proc tnsv(truth_vcf:string, pop_vcf:string, output_vcf:string="/dev/null", min_dist:int=100, max_tn_variants:int= -1) =
  var
    tvcf:VCF
    pvcf:VCF
  if not tvcf.open(truth_vcf):
    quit &"couldn't open truth_vcf: {truth_vcf}"
  if not pvcf.open(pop_vcf):
    quit &"couldn't open pop_vcf: {pop_vcf}"
  let tvcf_samples = tvcf.samples
  pvcf.set_samples(@["^"])

  #for sample in tvcf.samples:
  #  doAssert pvcf.header.hdr.bcf_hdr_add_sample(sample) >= 0
  discard pvcf.header.hdr.bcf_hdr_sync()

  var th = tvcf.header.hdr
  var ph = pvcf.header.hdr

  var chrs = initHashSet[string]()

  var mh = th.bcf_hdr_merge(ph)
  doAssert mh != nil, "error merging headers"

  # save variants when max_tn_variants is set.
  var tn_variants:seq[Variant]

  var ovcf:VCF
  if not ovcf.open(output_vcf, mode="w"):
    quit &"couldn't open output_vcf: {output_vcf}"
  ovcf.header = Header(hdr:mh)
  doAssert ovcf.write_header

  # read truth vcf into lapper
  var truth = newTable[string, seq[expanded_sv]]()
  for variant in tvcf:
    discard truth.hasKeyOrPut(stripChr(variant.CHROM), newSeq[expanded_sv]())
    var e = expanded_sv(start: max(0, variant.start.int - min_dist), stop: variant.stop.int + min_dist)
    truth[stripChr(variant.CHROM)].add(e)
    doAssert ovcf.write_variant(variant)
    chrs.incl($variant.CHROM)

  var tr = newTable[string, Lapper[expanded_sv]]()
  for chrom, evs in truth.mpairs:
    tr[chrom] = lapify(evs)

  # now iterate over pop and add any variant that doesn't overlap.
  var res = newSeq[expanded_sv]()
  var seen = initHashSet[string]()
  var overlapping = 0
  var new_variants = 0
  for v in pvcf:
    res.setLen(0)

    try:
      tr[stripChr(v.CHROM)].find(v.start.int, v.stop.int, res)
    except KeyError:
      continue

    if res.len != 0:
      overlapping += 1
      continue

    var key = v.key
    if key in seen: continue
    seen.incl(key)

    # new TN variant to add
    var vs = v.tostring().strip().split("\t")
    if $v.CHROM notin chrs:
      # try to match chroms between truth and pop set.
      var vc = $v.CHROM
      if vc.startsWith("chr") and stripChr(vc) in chrs:
        vs[0] = stripchr(vc)
      elif "chr" & vc in chrs:
        vs[0] = "chr" & vc

    if vs.len >= 10:
      vs = vs[0..<8]

    # HACK since I don't know how to add new fields to variant record
    if tvcf_samples.len > 0:
      vs.add("GT")
    for s in tvcf_samples:
      vs.add("0/0")
    var s = vs.join("\t")
    var ks = kstring_t(s:s, l:s.len.csize_t, m:s.len.csize_t)
    doAssert 0 == vcf_parse(ks.addr, mh, v.c)
    new_variants += 1

    if max_tn_variants > 0:
      tn_variants.add(v.copy())
    else:
      doAssert ovcf.write_variant(v)

  stderr.write_line "overlapping variants (not added):", overlapping

  if tn_variants.len > 0:
    shuffle(tn_variants)
    tn_variants = tn_variants[0..<min(tn_variants.len, max_tn_variants)]
    for v in tn_variants:
      doAssert ovcf.write_variant(v)
  else:
      stderr.write_line "new variants added:", new_variants

  ovcf.close


proc main() =

  var p = newParser("tnsv"):
    help("add population SVs as realistic true negatives (0/0 genotypes) to an existing truth-set")
    option("-m", "--min-dist", default="100", help="only include SVs from `pop_vcf` that are at least this far from the closes variant in `truth_vcf`")
    option("-n", "--max-tn-variants", default="-1", help="maximum number of TN variants to add (0 or less means add all from the truth-vcf)")
    #option("-i", "--include-types", default="CNV,INS,INV,DEL,DUP", help="only include SVs from `pop_vcf` that have one of these SVTYPES")
    option("-o", "--output-vcf", default="/dev/stdout", help="path for output truth vcf with new 0/0 variants added")
    arg("truth_vcf", help="VCF of SV truth-set")
    arg("pop_vcf", help="population VCF")
  try:
    var opts = p.parse()
    tnsv(opts.truth_vcf, opts.pop_vcf, output_vcf=opts.output_vcf, min_dist=parseInt(opts.min_dist), max_tn_variants=parseInt(opts.max_tn_variants))
  except UsageError as e:
    #stderr.write_line($e)
    stderr.write_line(p.help)
    stderr.write_line(getCurrentExceptionMsg())
    quit(1)

when isMainModule:
  main()
