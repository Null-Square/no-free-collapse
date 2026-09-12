# Paper and review package

The primary manuscript is **Complementary-subset dilations and a spectral hierarchy for hafnian energies**, in `matching_dilation.tex`. Build it twice with pdflatex from this directory; its only local input is `author_metadata.tex`. The generated PDF is an author-review candidate, not a submitted article.

## Frozen scientific scope

The paper proves the sharp complementary-subset operator dilation, the averaged spectral hierarchy for r>=2 and N>=4r-2, its hybrid with Roos, a sharp fixed-dimension density estimate, and all six-variable PSD equality cases. The real collision-free GBS result is a corollary, not a demonstrated experimental advantage.

The submission audit adds the missing odd-ambient-dimension justification for Roos, handles zero-energy cases explicitly, expands the equality proof, and strengthens the density ceiling for odd N. These changes are recorded in [the audit](../docs/SUBMISSION_AUDIT.md).

The complex zero-diagonal stability work remains a companion note and is not needed to prove the primary manuscript. The cube-normalized 1/216 problem remains open and is explicitly excluded from the submitted claim set.

## Author and submission decisions

No model or project label is assigned human authorship. The metadata input is intentionally marked pending. Confirm the human author list, order, affiliations, corresponding author, contribution statement, funding and competing interests. Read and approve the proofs and references, and approve the detailed AI-research and manuscript disclosure before changing the approval record.

See [submission/README.md](../submission/README.md). `python tools/submission_check.py` checks the technical package; `python tools/submission_check.py --submission` must remain nonzero until the recorded approvals are complete. A green CI is not human scientific sign-off.
