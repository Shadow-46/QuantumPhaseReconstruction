# Theory and Research Progress

## Research Objective

The goal of this research is to improve the classical reconstruction stage of Windowed Quantum Phase Estimation (WQPE) used inside the Modular Shor's Algorithm.

The objective is NOT to modify the quantum circuit itself.

Instead, the objective is to improve the classical reconstruction pipeline after quantum measurements have been obtained.

The long-term goal is to design a reconstruction framework that is more robust under ambiguity, limited resources, and noisy environments while maintaining shallow quantum circuits.

---

# Theoretical Background

## Standard Shor's Algorithm

I have studied Standard Shor's Algorithm in detail.

Topics understood include:

- Integer factorization problem
- Choosing N and a
- Coprime requirement
- Modular exponentiation
- Period finding
- Quantum Fourier Transform
- Inverse Quantum Fourier Transform
- Continued Fraction Expansion
- Recovering the order r
- Classical factor extraction

I understand the complete mathematical workflow from choosing a until recovering the factors.

---

## Quantum Phase Estimation

I understand:

- Eigenvalues
- Eigenstates
- Unitary operators
- Phase encoding
- Controlled-U operations
- Why phases are hidden before IQFT
- How IQFT converts phase information into measurable amplitudes
- Binary phase representation
- Phase precision
- Measurement probabilities
- Spectral leakage
- Dirichlet kernel intuition
- Relationship between phase estimation and period finding

---

## Windowed Quantum Phase Estimation

I understand why standard QPE is difficult for NISQ hardware.

I understand the motivation behind Windowed QPE:

- Smaller quantum circuits
- Lower circuit depth
- Reduced qubit requirements
- Multiple independent windows
- Classical reconstruction replacing one large IQFT

I understand:

- Window decomposition
- Window size
- Overlap bits
- Candidate generation
- Carry-aware stitching
- Reconstruction pipeline

---

## Modular Shor's Algorithm

I understand the complete modular algorithm proposed in the paper.

Including:

- Phase decomposition
- Independent window execution
- Candidate extraction
- Overlap consistency
- Carry computation
- Path stitching
- Beam-search-like reconstruction
- Continued Fraction recovery
- Final factor extraction

---

# Literature Review

I completed a comprehensive literature review.

The review covered:

- Standard QPE
- AWQPE
- Modular Shor
- Bayesian reconstruction
- Maximum Likelihood reconstruction
- Machine Learning based reconstruction
- Dynamic Programming
- Beam Search
- Graph Search
- Confidence-based decoding

The review also compared my original research idea against all existing work.

Major conclusions:

- Multi-candidate reconstruction already exists.
- Carry-aware stitching already exists.
- Beam-search-like pruning already exists.
- Static overlap already exists.
- Adaptive overlap has not been proposed.
- Continuous confidence-guided adaptive reconstruction appears to be novel.

---

# Repository Implementation

I implemented the complete published algorithm from scratch using Qiskit.

The repository includes:

- Quantum circuits
- IQFT
- QPE
- Modular multiplication
- Window generation
- Candidate extraction
- Carry-aware stitching
- Reconstruction pipeline
- Evaluation framework

The implementation reproduces the published algorithm.

---

# Verification

The implementation has been verified.

Each major module was individually tested.

The repository reproduces the published examples correctly.

The implementation is considered functionally correct.

---

# Experimental Campaign

A comprehensive experimental campaign has already been completed.

Experiments include:

- Parameter sensitivity
- Window size variation
- Overlap variation
- Beam width variation
- Candidate count variation
- Precision variation
- Noise injection
- Readout errors
- Depolarizing noise
- Phase damping
- Amplitude damping
- Thermal relaxation
- Scalability analysis
- Stress testing
- Adversarial testing
- Combined resource starvation
- Failure attribution

---

# Experimental Findings

The most important findings are:

1.

The published algorithm is robust under individual parameter variation.

2.

Failures mainly occur when multiple resource constraints interact.

3.

Carry validation accepts many incorrect candidate combinations.

4.

Ambiguous measurement regions significantly increase reconstruction uncertainty.

5.

Compilation complexity grows rapidly for larger modular multiplication circuits.

6.

The dominant algorithmic weakness appears during classical reconstruction under combined uncertainty rather than during quantum execution.

---

# Current Research Position

Implementation:
Completed

Verification:
Completed

Paper Reproduction:
Completed

Experimental Validation:
Completed

Failure Analysis:
Completed

The project is now entering the research design stage.

The next contribution should be motivated entirely by experimental evidence rather than intuition.

No new algorithm has been designed yet.

The goal is first to understand the failure landscape before proposing an improved reconstruction framework.

---

# Current Research Philosophy

The next contribution should:

- Preserve the published quantum circuit whenever possible.
- Modify only the classical reconstruction stage.
- Be supported by theoretical reasoning.
- Be supported by experimental evidence.
- Be statistically validated.
- Be compared directly against the published algorithm.

The objective is to produce a publishable improvement suitable for journals such as MethodsX followed by a higher-impact quantum algorithms journal.
