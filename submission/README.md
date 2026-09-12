# Submission package and remaining approvals

**Stage: technically prepared author-review candidate; not submitted.**

The next step is human proof/priority review and author approval, not further open-ended theorem hunting. The core result and its corollaries have complete written arguments. Independent human review cannot be replaced by another model-generated assurance.

## Proposed journal

The working first target is **Linear Algebra and its Applications**, chosen for the matrix-inequalities and finite-dimensional operator focus. This is a fit judgment, not a prediction of acceptance. Its official scope was checked, but the detailed current guide returned an access error in this pass. Consequently `venue_requirements_verified` remains false. Before upload, confirm its current source-file, abstract, keyword, highlights, declaration and supplementary-material requirements in the live submission portal.

The official Electronic Journal of Linear Algebra instructions were also checked as an alternative; they explicitly require its SIAM-derived LaTeX template, an abstract at most 250 words and 4-6 keywords. The supplied article-class candidate is not falsely labeled ELA-template compliant. No paid open-access choice has been made.

## Included materials

- `../paper/matching_dilation.tex` and `author_metadata.tex`: complete primary source, with unconfirmed author fields visibly marked.
- `cover_letter.md`: journal-addressed draft, with approval-dependent statements left unasserted.
- `highlights.txt`: four concise contribution statements.
- `ai_disclosure.md`: research and manuscript disclosure instructions reflecting the substantive AI use.
- `literature_review.md`: exact comparison points, sources examined and unresolved access/priority limitations.
- `review_record.json`: machine-readable author/declaration/review approvals, initially unset.
- `../tools/submission_check.py`: technical and strict submission preflights.

## Gates

Run `python tools/submission_check.py` to validate the technical inventory, LaTeX references and highlights. Run `python tools/submission_check.py --submission` to require author data, declarations, proof/priority sign-off, venue verification and exclusive-submission confirmation. Its nonzero exit while approvals are absent is intentional.

No human review, originality warranty, absence of conflicts, lack of funding, or author consent is fabricated. The author record and visible manuscript notices must be updated together only after actual approval. The final manuscript must accurately disclose AI's roles in theorem exploration, proof proposals, coding, literature triage and drafting, not merely proofreading.

## Final handoff sequence

Confirm human authors, order, affiliations and corresponding contact; arrange and record proof/priority review; finalize contributions, funding, interests and AI disclosure; check the journal's live requirements; replace the pending metadata and declarations; rebuild twice and visually check the final PDF; archive the exact approved source and commit; run the strict preflight; then obtain final authorization to submit. No journal submission or preprint upload has been performed.

Policy sources checked 12 September 2026: https://www.sciencedirect.com/journal/linear-algebra-and-its-applications/publish/guide-for-authors ; https://www.elsevier.com/about/policies-and-standards/generative-ai-policies-for-journals ; https://journals.uwyo.edu/index.php/ela/about/submissions .
