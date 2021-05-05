import hts/vcf
import random
import os
import strformat

proc main() =
  var ivcf:VCF
  if not ivcf.open(paramStr(1)): quit "couldn't open vcf"

  var tp = 0
  var tn = 0
  var fp = 0
  var fn = 0
  var neg = 0
  var pos = 0
  randomize()

  for v in ivcf:
    #if v.FILTER notin ["LongReadHomRef", "PASS", "."]: continue
    var method_calls_it_good = rand(1.0) < 0.96
    var x:seq[int32]
    var alts = v.format.genotypes(x).alts[0]
    if alts < 0: continue
    neg += int(alts == 0)
    pos += int(alts > 0)

    tp += int(method_calls_it_good and alts > 0)
    tn += int((not method_calls_it_good) and alts == 0)
    fn += int((not method_calls_it_good) and alts > 0)
    fp += int(method_calls_it_good and alts == 0)

  echo &"positive variants: {pos}"
  echo &"negative variants: {neg}"
  echo &"tp:{tp} fp:{fp} tn:{tn}"
  echo &"tp:{tp} fp:{fp} tn:{tn}"
  echo &"Recall (tp / (tp + fn)): { tp / (tp + fn):.2f}"
  echo &"Precision (tp / (tp + fp)): { tp / (tp + fp):.2f}"
  echo &"FDR (fp / (fp + tp)): {fp / (fp + tp):.2f}"

when isMainModule:
  main()
