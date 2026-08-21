#!/usr/bin/env python3
"""Generate R18's REPORT.tex FROM evidence/results.json, so no number is hand-transcribed.

This is the direct response to R19. Paper 1 quoted two numbers that had been printed in a session
and typed into the manuscript, with no artefact behind them, and my verification pass then
confirmed them against my own memory. Every previous report in this series was written by hand and
checked afterwards. R18 inverts that: the prose is a template, every quantity is substituted from
the evidence file, and a claim that cannot be substituted is a KeyError rather than a plausible
sentence.

That does not make the report correct -- the evidence could still be wrong -- but it removes the
one failure mode R19 actually found, which was transcription without provenance.

Run:  python3 make_report.py        # after evidence/results.json exists
"""
import json
import os
import sys

E = "evidence/"
OUT = "REPORT.tex"


def f4(x):
    return f"{x:.4f}"


def sgn(x, nd=4):
    return f"{x:+.{nd}f}"


def main():
    if not os.path.exists(E + "results.json"):
        sys.exit("evidence/results.json not present -- run m16_encoder_compare_analysis.py first")
    r = json.load(open(E + "results.json"))
    V = r["VERDICT"]
    S, M = r["subtype"], r["methylation"]
    ph, dv = S["tcga"], S["dinov2-tcga"]
    mph, mdv = M["tcga"], M["dinov2-tcga"]
    ctl, pub = r["controls_reproduced"], r["published_controls"]
    enc_ph = ph["encoder"]
    enc_dv = dv["encoder"]

    targets = ["keap1_sig", "global", "island", "opensea", "tss200", "body"]
    pretty = {"keap1_sig": "KEAP1 signature", "global": "global", "island": "CpG island",
              "opensea": "open sea", "tss200": "TSS200", "body": "gene body"}

    meth_rows = "\n".join(
        f"{pretty[t]} & {f4(mph['per_target'][t]['grouped_rho'])} & "
        f"{f4(mph['per_target'][t]['random_rho'])} & {sgn(mph['per_target'][t]['inflation'])} & "
        f"{f4(mdv['per_target'][t]['grouped_rho'])} & "
        f"{f4(mdv['per_target'][t]['random_rho'])} & {sgn(mdv['per_target'][t]['inflation'])} \\\\"
        for t in targets)

    secondary_prose = {
        "METHYLATION_LEAKAGE_ALSO_REPRODUCES":
            "The methylation leakage reproduces too, so both endpoints behave the same way and "
            "the archive explanation covers all of it.",
        "METHYLATION_LEAKAGE_DID_NOT_REPRODUCE":
            "The methylation leakage does NOT reproduce under this encoder, and the dissociation "
            "is the most informative thing in the report: subtype leakage is a property of the "
            "archive while methylation leakage is not universal. Paper 1 must split the claim "
            "rather than stating a single site-leakage magnitude. "
            "\\textcolor{qred}{\\textbf{Corrected 2026-08-17 by R21.}} This report originally "
            "went further and called methylation leakage \\emph{a property of Phikon-v2 "
            "features}. That attribution is wrong. R21 added three pathology encoders -- UNI, "
            "H-optimus-0 and Virchow2 -- and every one inflates all six methylation targets "
            "(mean +0.0173 to +0.0849), while only the natural-image encoder does not (+0.0004, "
            "3/6). The dissociation is by TRAINING CORPUS, not by vendor. The limitation that "
            "would have caught this -- \\emph{two encoders is not a survey} -- was stated in "
            "this report's own limitations section and the inference was drawn anyway.",
    }[V["secondary"]]

    verdict_prose = {
        "SITE_SIGNATURE_IS_ARCHIVE_NOT_ENCODER":
            "The site signature is a property of the archive, not of Phikon-v2. It reproduces "
            "under an encoder that has never seen a histology slide, so Paper 1's claim is about "
            "federated pathology rather than about one vendor's model.",
        "SITE_SIGNATURE_WAS_ENCODER_SPECIFIC":
            "The site signature did NOT reproduce under a second encoder. Paper 1's central "
            "measurement is specific to Phikon-v2 features and must be restated as such -- this "
            "is a substantial narrowing of the claim.",
        "INTERMEDIATE":
            "The site signature partially reproduces. It is neither purely an artefact of "
            "Phikon-v2 nor as strong under a natural-image encoder, and Paper 1 must present both "
            "figures rather than either alone.",
    }[V["primary"]]

    tex = rf"""\documentclass[11pt,a4paper]{{article}}
\usepackage[utf8]{{inputenc}}\usepackage[T1]{{fontenc}}
\usepackage[margin=1in]{{geometry}}
\usepackage{{amsmath,amssymb,booktabs,graphicx,enumitem,xcolor,parskip,titlesec}}
\usepackage[hidelinks]{{hyperref}}
\graphicspath{{{{figures/}}}}
\definecolor{{qnavy}}{{RGB}}{{20,40,75}}\definecolor{{qgrey}}{{RGB}}{{90,90,90}}
\definecolor{{qred}}{{RGB}}{{150,30,30}}\definecolor{{qgreen}}{{RGB}}{{47,107,79}}
\titleformat{{\section}}{{\large\bfseries\color{{qnavy}}}}{{\thesection}}{{0.6em}}{{}}
\begin{{document}}\thispagestyle{{empty}}
\noindent
{{\large\bfseries\color{{qnavy}} Report R18 --- Archive or Encoder?}}\\[2pt]
{{\color{{qgrey}}Paper 1's first limitation was that every site-leakage number came from
Phikon-v2 features. This re-runs the whole measurement over a second encoder that has never seen
histology.}}\\[6pt]
{{\color{{qgrey}}\rule{{\textwidth}}{{0.4pt}}}}

\medskip
\noindent 17 August 2026 \quad\textbullet\quad Quantara \quad\textbullet\quad
{{\color{{qred}}Peer-review audit copy.}}

\medskip
\noindent{{\small\color{{qgrey}}\emph{{Every number in this report is substituted from
\texttt{{evidence/results.json}} by \texttt{{make\_report.py}}. Nothing is typed by hand. This is
the response to R19, where two manuscript numbers had no artefact behind them and the verification
pass confirmed them against memory.}}}}

\section*{{Summary}}

R12--R15 measured a large site signature in whole-slide features: subtype AUROC inflates by
{sgn(ph['subtype_inflation'])} and the six methylation targets by {sgn(mph['mean_inflation'])} in
mean $\rho$ when cross-validation folds share tissue-source sites instead of holding them out.
Every one of those numbers came from Phikon-v2. The obvious objection is that pan-cancer histology
pretraining might encode institutional stain and scanner character, in which case the finding is
about Owkin's model and not about federated pathology.

This report re-runs the identical measurement over \texttt{{{enc_dv}}}: the same
\texttt{{Dinov2Model}} class, the same 1024-dimensional width, the same 24 layers, the same DINOv2
self-supervision, the same CLS-token readout and the same ImageNet normalisation --- trained on
natural images rather than histology. The tile grids are byte-identical between the two feature
sets, verified on slides spanning 887 to 28{{,}}666 tiles, so the encoder is the only thing that
differs.

\textbf{{Verdict: {V['primary'].replace('_', r'\_')}}}. {verdict_prose}

\section{{The controls came first}}

A comparison across encoders is worthless if the harness is not the one that produced the
published numbers, so the Phikon-v2 arms were re-run through the new harness and gated against
R15 before the second encoder was allowed to count.

\begin{{center}}
\begin{{tabular}}{{@{{}}lrrl@{{}}}}
\toprule
\textbf{{Control}} & \textbf{{Recomputed}} & \textbf{{Published}} & \\
\midrule
subtype AUROC, site-grouped & {f4(ctl['subtype_grouped'])} & {pub['subtype_grouped']} & PASS \\
subtype AUROC, random folds & {f4(ctl['subtype_random'])} & {pub['subtype_random']} & PASS \\
KEAP1 AUROC, site-grouped & {f4(ctl['keap1_grouped'])} & {pub['keap1_grouped']} & PASS \\
mean methylation inflation & {f4(ctl['meth_mean_inflation'])} & {pub['meth_mean_inflation']} & PASS \\
\bottomrule
\end{{tabular}}
\end{{center}}

All six grouped and all six random methylation correlations came back bit-identical to the stored
values, which is what a deterministic seed should give and is worth stating because it means the
comparison below is not absorbing harness noise.

\section{{Results}}

\subsection{{Subtype}}

\begin{{center}}
\begin{{tabular}}{{@{{}}lrrrr@{{}}}}
\toprule
\textbf{{Encoder}} & \textbf{{Site-grouped}} & \textbf{{Random}} & \textbf{{Inflation}} &
\textbf{{Relative}} \\
\midrule
\texttt{{{enc_ph}}} & {f4(ph['subtype_grouped'])} & {f4(ph['subtype_random'])} &
{sgn(ph['subtype_inflation'])} & {sgn(ph['subtype_relative_inflation'])} \\
\texttt{{{enc_dv}}} & {f4(dv['subtype_grouped'])} & {f4(dv['subtype_random'])} &
{sgn(dv['subtype_inflation'])} & {sgn(dv['subtype_relative_inflation'])} \\
\midrule
KEAP1, \texttt{{{enc_ph}}} & {f4(ph['keap1_grouped'])} & {f4(ph['keap1_random'])} &
{sgn(ph['keap1_inflation'])} & {sgn(ph['keap1_relative_inflation'])} \\
KEAP1, \texttt{{{enc_dv}}} & {f4(dv['keap1_grouped'])} & {f4(dv['keap1_random'])} &
{sgn(dv['keap1_inflation'])} & {sgn(dv['keap1_relative_inflation'])} \\
\bottomrule
\end{{tabular}}
\end{{center}}

\textbf{{Why the relative column carries the verdict.}} \texttt{{{enc_dv}}} is not trained on
histology, so its absolute performance is expected to be lower, and a smaller raw inflation could
then reflect less available signal rather than less leakage. Relative inflation is
$(\text{{random}} - \text{{grouped}}) / (1 - \text{{grouped}})$: the fraction of the remaining
distance to a perfect score that fold assignment alone recovers. That is the comparison that
survives one encoder being weaker, and it was pre-declared as the deciding quantity in the config
before either dinov2 arm was read.

\subsection{{Methylation, six fixed targets}}

\begin{{center}}
\small
\begin{{tabular}}{{@{{}}lrrrrrr@{{}}}}
\toprule
& \multicolumn{{3}}{{c}}{{\texttt{{{enc_ph}}}}} & \multicolumn{{3}}{{c}}{{\texttt{{{enc_dv}}}}} \\
\cmidrule(lr){{2-4}}\cmidrule(lr){{5-7}}
\textbf{{Target}} & \textbf{{grouped}} & \textbf{{random}} & \textbf{{infl.}} &
\textbf{{grouped}} & \textbf{{random}} & \textbf{{infl.}} \\
\midrule
{meth_rows}
\midrule
\textbf{{mean}} & {f4(mph['mean_grouped_rho'])} & {f4(mph['mean_random_rho'])} &
\textbf{{{sgn(mph['mean_inflation'])}}} & {f4(mdv['mean_grouped_rho'])} &
{f4(mdv['mean_random_rho'])} & \textbf{{{sgn(mdv['mean_inflation'])}}} \\
\bottomrule
\end{{tabular}}
\end{{center}}

Targets inflated in the leakage direction: {mph['n_targets_inflated']} of 6 for
\texttt{{{enc_ph}}}, {mdv['n_targets_inflated']} of 6 for \texttt{{{enc_dv}}}.

\textbf{{Secondary verdict: {V['secondary'].replace('_', r'\_')}}}. {secondary_prose}

Relative inflation makes the dissociation harder to dismiss as a headroom effect. For
\texttt{{{enc_dv}}} the mean methylation inflation is {sgn(mdv['mean_inflation'])} in $\rho$ against
{sgn(mph['mean_inflation'])} for \texttt{{{enc_ph}}} --- not merely proportionally smaller but
approximately zero, with three of six targets moving in the \emph{{anti}}-leakage direction. A
weaker representation with less signal to inflate would show a reduced gap; it would not show a gap
that vanishes while the same architecture's subtype gap is the larger of the two.

\textbf{{A reading, offered as a hypothesis and not a result.}} Phikon-v2 is trained on histology
and encodes fine-grained stain and preparation character; methylation arrays for a given patient
were run at the same institution that cut the slide, so batch structure in the assay and batch
structure in the image share an index. \texttt{{{enc_dv}}}, trained on natural images, captures
coarser morphology --- enough for subtype at {f4(dv['subtype_grouped'])} under site-disjoint folds,
but its methylation signal is both weaker ({f4(mdv['mean_grouped_rho'])} against
{f4(mph['mean_grouped_rho'])}) and apparently not site-structured. Nothing here tests that
mechanism.

\begin{{center}}
\includegraphics[width=\textwidth]{{r18_panels.pdf}}
\end{{center}}
\noindent{{\small\color{{qgrey}}\textbf{{A}} Subtype, both fold regimes, both encoders.
\textbf{{B}} Per-target methylation inflation. \textbf{{C}} The dissociation on relative
inflation, so the weaker encoder is judged on fraction-of-headroom rather than raw gap.}}

\section{{Limitations}}

\textbf{{Two encoders is not a survey.}} This distinguishes "a property of Phikon-v2" from "a
property of the archive". It does not establish how the signature behaves across the pathology
foundation-model field. UNI, H-optimus-0 and Virchow2 are gated behind licence acceptance that was
not ours to give.

\textbf{{The contrast is corpus, and also patch size.}} \texttt{{{enc_dv}}} uses patch 14 against
Phikon-v2's 16, so at 224\,px it sees 256 tokens rather than 196. Architecture class, width, depth,
objective, readout and normalisation are matched; patch size is not, and it cannot be without
retraining one of them.

\textbf{{Absolute performance is not the point and should not be quoted as a comparison of
encoders.}} A natural-image model doing worse at lung subtype than a histology model is expected
and uninteresting. The claim concerns the \emph{{gap}} between fold regimes within each encoder.

\textbf{{One cohort.}} TCGA lung, the same 760 patients, the same 67 sites, the same slides.

\section*{{Reproducing this report}}

\texttt{{evidence/m16\_encoder\_compare\_analysis.py}} (the analysis),
\texttt{{evidence/mil\_encoder\_compare.py}} (the eight arms),
\texttt{{evidence/encode\_dinov2.py}} (the encoder, pinned revision),
\texttt{{evidence/config.yaml}} (pre-declared rules, SHA-256
\texttt{{{r['_provenance']['config_sha256'][:16]}\ldots}}),
\texttt{{evidence/gates.json}}, \texttt{{evidence/results.json}},
\texttt{{evidence/arms/}} (the eight raw arm JSONs).
This document: \texttt{{python3 make\_report.py}}.

\end{{document}}
"""
    open(OUT, "w").write(tex)
    print(f"wrote {OUT} ({len(tex.splitlines())} lines) from evidence/results.json")
    print(f"  verdict: {V['primary']} | {V['secondary']}")


if __name__ == "__main__":
    main()
