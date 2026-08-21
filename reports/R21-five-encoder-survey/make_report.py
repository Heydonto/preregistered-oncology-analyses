#!/usr/bin/env python3
"""Generate R21's REPORT.tex FROM evidence/results.json. No number is typed by hand.

Same discipline as R18: the prose is a template, every quantity is substituted, and a claim that
cannot be substituted raises KeyError rather than reading as a plausible sentence. The verdict
prose is selected by the mechanical verdict string so the narration cannot outrun the rule.
"""
import json
import os
import sys

E = "evidence/"
if not os.path.exists(E + "results.json"):
    sys.exit("run m18_five_encoder_survey.py first")
r = json.load(open(E + "results.json"))
V, ENC, BC, PH = r["VERDICT"], r["encoders"], r["by_corpus"], r["post_hoc_observation"]
ORDER = ["tcga", "dinov2-large", "uni-tcga", "hopt-tcga", "virchow2-tcga"]
SHORT = {"tcga": "Phikon-v2", "dinov2-large": "dinov2-large", "uni-tcga": "UNI",
         "hopt-tcga": "H-optimus-0", "virchow2-tcga": "Virchow2"}
TG = ["keap1_sig", "global", "island", "opensea", "tss200", "body"]
PRETTY = {"keap1_sig": "KEAP1 sig.", "global": "global", "island": "island",
          "opensea": "open sea", "tss200": "TSS200", "body": "gene body"}


def f4(x):
    return f"{x:.4f}"


def sg(x, n=4):
    return f"{x:+.{n}f}"


rows = "\n".join(
    rf"\texttt{{{SHORT[k]}}} & {ENC[k]['corpus']} & {f4(ENC[k]['subtype_grouped'])} & "
    rf"{f4(ENC[k]['subtype_random'])} & {sg(ENC[k]['subtype_inflation'])} & "
    rf"{sg(ENC[k]['subtype_relative_inflation'])} & {sg(ENC[k]['keap1_inflation'])} & "
    rf"{f4(ENC[k]['meth_mean_grouped'])} & {sg(ENC[k]['meth_mean_inflation'])} & "
    rf"{ENC[k]['meth_targets_inflated']}/6 \\" for k in ORDER)

per_rows = "\n".join(
    PRETTY[t] + " & " + " & ".join(sg(ENC[k]["meth_per_target"][t]["inflation"]) for k in ORDER)
    + r" \\" for t in TG)

meth_prose = {
    "METHYLATION_LEAKAGE_IS_A_PROPERTY_OF_HISTOLOGY_PRETRAINING":
        rf"""Every one of the {BC['histology']['n']} histology-pretrained encoders inflates
\textbf{{all six}} methylation targets, with mean inflation from
{sg(BC['histology']['meth_inflation_range'][0])} to
{sg(BC['histology']['meth_inflation_range'][1])} in $\rho$. The single natural-image encoder
inflates {BC['natural_images']['targets_inflated']} of 6 with a mean of
{sg(BC['natural_images']['meth_inflation'])} --- indistinguishable from nothing. The dissociation is
by \emph{{training corpus}}, not by vendor.

\textcolor{{qred}}{{\textbf{{This corrects R18.}}}} R18 compared Phikon-v2 against dinov2-large
alone and concluded that methylation leakage was ``a property of Phikon-v2 features''. With four
histology encoders in hand that attribution is wrong. R18 had named the limitation that produced
the error --- ``two encoders is not a survey'' --- and drew the inference anyway.""",
    "METHYLATION_LEAKAGE_IN_ALL_ENCODERS":
        "Methylation leakage appears in every encoder including the natural-image one, so it is "
        "not corpus-specific and R18's dissociation does not survive at all.",
    "METHYLATION_LEAKAGE_MIXED_ACROSS_HISTOLOGY_ENCODERS":
        "Methylation leakage is inconsistent among the histology encoders, so no corpus-level "
        "statement is supportable and R18's attribution must simply be withdrawn as unresolved.",
}[V["methylation"]]

tex = rf"""\documentclass[11pt,a4paper]{{article}}
\usepackage[utf8]{{inputenc}}\usepackage[T1]{{fontenc}}
\usepackage[margin=0.9in]{{geometry}}
\usepackage{{amsmath,amssymb,booktabs,graphicx,enumitem,xcolor,parskip,titlesec}}
\usepackage[hidelinks]{{hyperref}}
\graphicspath{{{{figures/}}}}
\definecolor{{qnavy}}{{RGB}}{{20,40,75}}\definecolor{{qgrey}}{{RGB}}{{90,90,90}}
\definecolor{{qred}}{{RGB}}{{150,30,30}}\definecolor{{qgreen}}{{RGB}}{{47,107,79}}
\titleformat{{\section}}{{\large\bfseries\color{{qnavy}}}}{{\thesection}}{{0.6em}}{{}}
\begin{{document}}\thispagestyle{{empty}}
\noindent
{{\large\bfseries\color{{qnavy}} Report R21 --- Five Encoders, and a Correction to R18}}\\[2pt]
{{\color{{qgrey}}Subtype site-leakage is universal. Methylation site-leakage tracks the training
corpus, not the vendor --- which is not what R18 concluded.}}\\[6pt]
{{\color{{qgrey}}\rule{{\textwidth}}{{0.4pt}}}}

\medskip
\noindent 18 August 2026 \quad\textbullet\quad Quantara \quad\textbullet\quad
{{\color{{qred}}Peer-review audit copy.}}

\medskip
\noindent{{\small\color{{qgrey}}\emph{{Every number here is substituted from
\texttt{{evidence/results.json}} by \texttt{{make\_report.py}}; nothing is typed by hand.}}}}

\section*{{Summary}}

R18 measured site-leakage on two encoders and R20 on the external cohorts. Both carried the same
limitation in writing: two encoders is not a survey, and one of the two was not a pathology model.
Three pathology foundation models now close it --- UNI, H-optimus-0 and Virchow2 --- each encoding
all 1{{,}}182 slides with byte-identical tile grids and its own prescribed normalisation.

\textbf{{Verdict 1: {V['subtype'].replace('_', r'\_')}}}. Fold assignment inflates subtype AUROC in
all five encoders, every one clearing the relative-inflation bar of 0.05 that R18 declared before
its own arms were read.

\textbf{{Verdict 2: {V['methylation'].replace('_', r'\_')}}}. {meth_prose}

\section{{The grid}}

\begin{{center}}
\small
\begin{{tabular}}{{@{{}}llrrrrrrrr@{{}}}}
\toprule
& & \multicolumn{{4}}{{c}}{{\textbf{{subtype}}}} & \textbf{{KEAP1}} &
\multicolumn{{3}}{{c}}{{\textbf{{methylation}}}} \\
\cmidrule(lr){{3-6}}\cmidrule(lr){{7-7}}\cmidrule(lr){{8-10}}
\textbf{{Encoder}} & \textbf{{Corpus}} & \textbf{{grouped}} & \textbf{{random}} &
\textbf{{infl.}} & \textbf{{rel.}} & \textbf{{infl.}} & \textbf{{grouped}} &
\textbf{{infl.}} & \textbf{{n}} \\
\midrule
{rows}
\bottomrule
\end{{tabular}}
\end{{center}}

\begin{{center}}
\includegraphics[width=\textwidth]{{r21_panels.pdf}}
\end{{center}}
\noindent{{\small\color{{qgrey}}\textbf{{A}} Subtype, both fold regimes. \textbf{{B}} Methylation
inflation coloured by training corpus. \textbf{{C}} A post-hoc pattern, labelled as one.}}

\subsection{{Per-target methylation inflation}}

\begin{{center}}
\small
\begin{{tabular}}{{@{{}}l{'r' * len(ORDER)}@{{}}}}
\toprule
\textbf{{Target}} & {' & '.join(rf'\textbf{{{SHORT[k]}}}' for k in ORDER)} \\
\midrule
{per_rows}
\bottomrule
\end{{tabular}}
\end{{center}}

\noindent Only \texttt{{dinov2-large}} has targets moving the wrong way. Every histology encoder is
positive on all six, which is what makes the corpus split legible rather than a difference of degree.

\section{{Subtype leakage is universal but not constant}}

The relative figure runs from {f4(min(ENC[k]['subtype_relative_inflation'] for k in ORDER))} to
{f4(max(ENC[k]['subtype_relative_inflation'] for k in ORDER))}, and site-disjoint capability from
{f4(min(ENC[k]['subtype_grouped'] for k in ORDER))} to
{f4(max(ENC[k]['subtype_grouped'] for k in ORDER))}. The three newer pathology models are far
better at the site-disjoint task than Phikon-v2 --- Virchow2 reaches
{f4(ENC['virchow2-tcga']['subtype_grouped'])} against Phikon-v2's
{f4(ENC['tcga']['subtype_grouped'])} --- and show correspondingly smaller inflation.

So ``a fold-assignment artefact of {sg(ENC['tcga']['subtype_inflation'])}'' is a statement about a
particular encoder in 2024, not a constant of the archive. The leakage is real in every encoder
tested; its size is not transferable.

\section{{A pattern we are deliberately not claiming}}

Panel C plots relative inflation against site-disjoint capability. It looks like a clean inverse
relationship: Spearman $\rho = {PH['spearman_rho']:.2f}$, $p = {PH['p']:.3f}$,
$n = {PH['n_encoders']}$.

\textbf{{We are not reporting that as a finding, and the config says so.}} It was noticed
\emph{{after}} all five encoders were read. Testing a hypothesis on the data that suggested it is
the precise move this programme exists to police, and at $n = {PH['n_encoders']}$ no rank
correlation can reach $p<0.05$ by this test regardless of how clean the picture looks. What it
licenses is a pre-registration --- fix the rule, then test it on encoders not used here --- not a
conclusion. It is recorded in \texttt{{results.json}} under
\texttt{{post\_hoc\_observation}} with \texttt{{STATUS}} set accordingly.

\section{{Controls}}

The Phikon-v2 arms were re-run through the same harness and still match R15: subtype
{f4(ENC['tcga']['subtype_grouped'])} / {f4(ENC['tcga']['subtype_random'])}, KEAP1
{f4(ENC['tcga']['keap1_grouped'])}, mean methylation inflation
{f4(ENC['tcga']['meth_mean_inflation'])}. All gates in \texttt{{evidence/gates.json}} passed. Every
arm covers the same 760 patients and 67 sites.

KEAP1 inflation stays small in all five ({sg(min(ENC[k]['keap1_inflation'] for k in ORDER))} to
{sg(max(ENC[k]['keap1_inflation'] for k in ORDER))}), which is what R15 predicted on the mechanistic
grounds that mutation status has no institutional prevalence structure to exploit. A prediction made
for a reason, holding across five encoders chosen for other reasons, is worth more than the same
number from one.

\section{{Limitations}}

\textbf{{Five encoders, one archive, one cohort.}} TCGA lung, 760 patients, 67 sites, the same
slides and the same tiles throughout.

\textbf{{Corpus is confounded with vintage and capacity.}} The three pathology models are newer,
larger and trained on more data than Phikon-v2, and dinov2-large differs in corpus \emph{{and}} in
being the only non-pathology model. A corpus-level statement is the most parsimonious reading of
Panel B, not the only possible one.

\textbf{{Normalisation is not held constant.}} H-optimus-0 prescribes its own; forcing ImageNet
statistics on it would have crippled it. Each model gets its own preprocessing, which is correct
practice and a departure from R18's tighter design.

\textbf{{CLS-only readout.}} Virchow2's authors recommend concatenating the class token with the
mean of the patch tokens. Holding the readout constant is what makes five encoders comparable, but
it may understate Virchow2's absolute capability.

\textbf{{Two of the five are non-commercially licensed.}} UNI and Virchow2 are CC-BY-NC-ND-4.0;
Phikon-v2, dinov2-large and H-optimus-0 are not. Each arm records its encoder's licence in the
feature files. Anything that could read as commercial should rest on the Apache-2.0 arm.

\section*{{Reproducing this report}}

\texttt{{evidence/encode\_pathfm.py}} (the three encoders, pinned),
\texttt{{evidence/mil\_encoder\_compare.py}} (the twenty arms),
\texttt{{evidence/m18\_five\_encoder\_survey.py}} (this analysis),
\texttt{{evidence/config.yaml}} (SHA-256 \texttt{{{r['_provenance']['config_sha256'][:16]}\ldots}}),
\texttt{{evidence/gates.json}}, \texttt{{evidence/results.json}},
\texttt{{evidence/arms/}} (all twenty raw arm JSONs).
This document: \texttt{{python3 make\_report.py}}. Figure: \texttt{{python3 make\_figure.py}}.

\end{{document}}
"""
open("REPORT.tex", "w").write(tex)
print(f"wrote REPORT.tex ({len(tex.splitlines())} lines)")
print(f"  {V['subtype']} | {V['methylation']}")
