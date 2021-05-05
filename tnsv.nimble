# Package

version       = "0.0.1"
author        = "Brent Pedersen"
description   = "add true-negatives to an SV truth-set"
license       = "MIT"

# Dependencies
requires "nim >= 0.19.2", "hts >= 0.3.15", "lapper", "argparse == 0.10.1"

bin = @["tnsv"]
