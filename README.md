This software adds true-negative SVs to a truth-set. The true-negatives should be drawn from a real set.

In short you can use like this:
```
wget https://github.com/brentp/tnsv/releases/download/v0.0.1/tnsv
chmod +x ./tnsv
./tnsv sv-truth-set.vcf.gz \
       population-sv-calls.vcf.gz \
    | bcftools sort -O z -o HG002_SVS.with-gnomad-TN.vcf.gz
```

And many **non-overlapping** hom-ref (genotype 0/0) calls will be added to the truth-set (in this case, the output is `HG002_SVS.with-gnomad-TN.vcf.gz`)
These calls will be in realistic locations (compared to random locations).

See [below](#true-negative-sets) for links to some possible population call-sets.

# Problem

This is a basic, known data-science problem, but it can still catch even seasoned analysts and it's especially
easy to hit this problem in genomics.

Given a truth-set like the [Genome in a Bottle SV set](https://www.nature.com/articles/s41587-020-0538-8)
we can evaluate methods for SV detection and filtering.

In `examples/svfilter.nim` I have built a sophisticated classifier that will evaluate a set of SVs
and **randomly retain (classify as valid) 96% of variants**.

This method is able to achieve:
+ 96% recall
+ 82% precision

on the actual [HG002 SV truth set](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz). 

# How!!!?

**recall** is (true-positives / (true-positives + false-negatives)), since we classify 96% of variants
as true, then our recall must be 96%.


**precision** is (true-positives / (true-positives + false-positives)). So why is precision so high? 
We have: `10757` negative variants and `48606` positive variants. Since the classifier is randomly choosing support
or not, then we can expect the precision is: `48606 / (48606 + 10757)` which gives us the **82%**.


In short, because there are so few negatives, a random classifier (with a bias toward the positives) will appear to have decent or even
great performance.

# Mitigation

One way to make it harder to miss this problem is to **add many true-negative** variants.
Instead of doing this randomly, we add true variants from a given population or database set.

Only population variants that are not within `dist` bases (default 100) of a variant in the truth-set are added.
This ensures that the added variants are realistic (not random locations).

For example, to add [gnomad-sv calls](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd166.GRCh37.variant_call.vcf.gz) to the Genome in a bottle truth set, use:
```
# tnsv $truthset $populationset > $augmented_truth_set
tnsv HG002_SVs_Tier1_v0.6.vcf.gz \
           nstd166.GRCh37.variant_call.vcf.gz \
	     | bcftools sort -O z -o HG002_SVS.with-gnomad-TN.vcf.gz
```

Now we can re-try our random classifer with the following results:

+ 96% Recall (which must be the case)
+ 15% Precision (down from 82% on the original truth-set)

With this, we have a such a low precision that we should note that something is wrong.

# True Negative sets:

+ [gnomad-svs GRCh37](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd166.GRCh37.variant_call.vcf.gz)
+ [gnomad-svs GRCh38](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd166.GRCh38.variant_call.vcf.gz)
+ [clinvar svs GRCh37](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd102.GRCh37.variant_call.vcf.gz)
+ [clinvar svs GRCh38](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd102.GRCh38.variant_call.vcf.gz)

`tnsv` will do simple re-mapping of chromosomes to match the truth-set so that e.g. '22' in the population set can become 'chr22'
if 'chr22' is present in the truth-set.

