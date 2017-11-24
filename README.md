<!-- Added by Tom:  mathjax -->
<script type="text/javascript" async
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML">
</script>


Phase retrieval is the recovery of a signal from only the magnitudes, and not the phases, of complex-valued linear measurements.  Phase retrieval problems arise in many different applications,
particularly in crystallography and microscopy.
Mathematically, phase retrieval recovers a complex valued signal \\(x\in \mathbb{C}^n\\) from \\(m\\) measurements of the form

\\begin{align} \label{eq:original}
b_i = | a_i^T x|, \quad i = 1,2,\ldots,m,
\\end{align}

where \\(\\{a_i\\}\\) are complex-valued measurement vectors.

 PhasePack is a collection of sub-routines for solving classical phase retrieval problems.  PhasePack contains implementations of both classical and contemporary phase retrieval routines.  It also contains scripts that apply these signal recovery
 routines to datasets, and compares the results of different methods in terms of speed and quality.  PhasePack can do comparisons not only with synthetic datasets, but also with several publicly available real-world datasets.

 The routines implemented in PhasePack share a simple common interface, making it easy to switch between
 solvers.  The interface also gives users the option to control the initialization and runtime of different algorithms, making it easy to "plug and play" different routines and options.
