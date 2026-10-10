"""dmpkCTGcounter: CTG repeat units of the two DMPK alleles from ONT reads.

The package takes one BAM or CRAM (GRCh38, ONT long reads) and reports the
(CTG)n repeat count of the DMPK 3' UTR for the two alleles as median and
range. Coordinates in the public interface are 1-based and closed (SAM and
VCF convention); internal sequence slicing is 0-based and half-open, and the
two are converted where the reference is read (``locus.verify_reference``)
and where records are fetched (``pipeline.count_ctg``).
"""

import os

__version__ = "0.1.0"

# htslib looks up a CRAM reference by MD5 in REF_PATH and, by default, in a
# public registry over the network. Point both variables at paths that do not
# exist so that no lookup, cache, or download happens. This does not disable
# the local @SQ UR path that a CRAM header may carry: if that file still
# exists, htslib uses it. A reference given with --reference is tried first.
os.environ["REF_PATH"] = "/nonexistent/dmpkctgcounter-no-ref-path"
os.environ["REF_CACHE"] = "/nonexistent/dmpkctgcounter-no-ref-cache"
