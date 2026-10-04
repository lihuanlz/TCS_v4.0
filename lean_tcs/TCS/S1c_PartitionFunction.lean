/-
=====================================================================
S1c: combinatorial core of the partition function
=====================================================================
Corresponds to SI S1c: the statistical-mechanics setup behind the master
equation. The physical correspondence (capture sites as a canonical
ensemble, Boltzmann weights) is a modelling assumption; its combinatorial
content is pure mathematics and is formalized here:

  - partitionFunction_factorizes: the grand partition function of Ω
    independent two-state sites, Z = Σ_k C(Ω,k)·x^k, equals (1+x)^Ω
    (binomial theorem, mathlib's Commute.add_pow)
  - occupancyNumerator_eq: the first-moment sum satisfies
    Σ_k k·C(Ω,k)·x^k = Ω·x·(1+x)^(Ω−1)
  - mean_occupancy: mean occupancy p = x/(1+x) (i.e. (x/Ω)·∂lnZ/∂x
    evaluated algebraically through the two sum identities)
  - occupancy_pmf: the occupied-site count is Binomial(Ω, p) at PMF
    level: C(Ω,k)·x^k/(1+x)^Ω = C(Ω,k)·p^k·(1−p)^(Ω−k)
  - activity_from_occupancy: inversion x = p/(1−p), i.e. the Langmuir
    relation (Axiom S2e.1), which is the κ→∞ limit of the master
    equation ξ = p/(1−p) + p/κ (Eq. S1c.12)

Lean 4 + Mathlib, toolchain leanprover/lean4:v4.15.0
-/

import Mathlib

open Finset

namespace TCS

/-- Grand partition function of Ω independent two-state capture sites at
activity x: configurations with k occupied sites carry degeneracy
C(Ω,k) and Boltzmann weight x^k. -/
noncomputable def partitionFunction (Ω : ℕ) (x : ℝ) : ℝ :=
  ∑ k ∈ Finset.range (Ω + 1), (Ω.choose k : ℝ) * x ^ k

/-- Occupancy numerator: the first moment of the configuration sum,
Σ_k k·C(Ω,k)·x^k. -/
noncomputable def occupancyNumerator (Ω : ℕ) (x : ℝ) : ℝ :=
  ∑ k ∈ Finset.range (Ω + 1), (k : ℝ) * (Ω.choose k : ℝ) * x ^ k

/-- Factorization of the partition function (binomial theorem):
Z = Σ_k C(Ω,k)·x^k = (1+x)^Ω. -/
theorem partitionFunction_factorizes (Ω : ℕ) (x : ℝ) :
    partitionFunction Ω x = (1 + x) ^ Ω := by
  unfold partitionFunction
  rw [add_comm (1 : ℝ) x, Commute.add_pow (Commute.one_right x) Ω]
  apply Finset.sum_congr rfl
  intro k _
  simp [one_pow, mul_comm, mul_left_comm]

/-- First-moment sum: Σ_k k·C(Ω,k)·x^k = Ω·x·(1+x)^(Ω−1).
Uses (k+1)·C(n+1,k+1) = (n+1)·C(n,k) (Nat.succ_mul_choose_eq) and the
factorization of the partition function at Ω−1. -/
theorem occupancyNumerator_eq (Ω : ℕ) (x : ℝ) :
    occupancyNumerator Ω x = (Ω : ℝ) * x * (1 + x) ^ (Ω - 1) := by
  cases Ω with
  | zero => simp [occupancyNumerator]
  | succ n =>
    have hterm : ∀ k : ℕ,
        ((k + 1 : ℕ) : ℝ) * (((n + 1).choose (k + 1) : ℕ) : ℝ) * x ^ (k + 1)
          = (((n + 1 : ℕ)) : ℝ) * x * ((((n.choose k : ℕ)) : ℝ) * x ^ k) := by
      intro k
      have h := congrArg (fun m : ℕ => (m : ℝ)) (Nat.succ_mul_choose_eq n k)
      simp only [Nat.cast_mul, Nat.succ_eq_add_one, Nat.cast_add,
        Nat.cast_one] at h
      have hcastk : ((k + 1 : ℕ) : ℝ) = (k : ℝ) + 1 := by simp
      have hcastn : ((n + 1 : ℕ) : ℝ) = (n : ℝ) + 1 := by simp
      rw [hcastk, hcastn, pow_succ]
      have hR : ((k : ℝ) + 1) * (((n + 1).choose (k + 1) : ℕ) : ℝ)
          = ((n : ℝ) + 1) * ((n.choose k : ℕ) : ℝ) := by linarith [h]
      rw [hR]
      ring
    have hsum : (∑ k ∈ Finset.range (n + 1),
          ((k + 1 : ℕ) : ℝ) * (((n + 1).choose (k + 1) : ℕ) : ℝ) * x ^ (k + 1))
        = (((n + 1 : ℕ)) : ℝ) * x
            * (∑ k ∈ Finset.range (n + 1), (((n.choose k : ℕ)) : ℝ) * x ^ k) := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro k _
      exact hterm k
    have hZ : (∑ k ∈ Finset.range (n + 1), (((n.choose k : ℕ)) : ℝ) * x ^ k)
        = (1 + x) ^ n :=
      partitionFunction_factorizes n x
    simp only [occupancyNumerator]
    rw [Finset.sum_range_succ']
    simp only [Nat.cast_zero, zero_mul, add_zero]
    rw [show n + 1 - 1 = n from Nat.add_sub_cancel n 1]
    show (∑ k ∈ Finset.range (n + 1),
        ((k + 1 : ℕ) : ℝ) * (((n + 1).choose (k + 1) : ℕ) : ℝ) * x ^ (k + 1))
        = (((n + 1 : ℕ)) : ℝ) * x * (1 + x) ^ n
    rw [hsum, hZ]

/-- Mean occupancy of Ω independent two-state sites: p = x/(1+x).
This is (x/Ω)·∂lnZ/∂x evaluated algebraically via the two sum
identities. -/
theorem mean_occupancy {Ω : ℕ} (hΩ : Ω ≠ 0) {x : ℝ} (hx : (1 : ℝ) + x ≠ 0) :
    occupancyNumerator Ω x / ((Ω : ℝ) * partitionFunction Ω x) = x / (1 + x) := by
  obtain ⟨m, rfl⟩ := Nat.exists_eq_succ_of_ne_zero hΩ
  simp only [Nat.succ_eq_add_one]
  rw [occupancyNumerator_eq, partitionFunction_factorizes,
    show m + 1 - 1 = m from Nat.add_sub_cancel m 1, pow_succ]
  have hm : (((m + 1 : ℕ)) : ℝ) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero m
  have hZ : (1 + x) ^ m ≠ 0 := pow_ne_zero m hx
  field_simp
  ring

/-- Inversion of the occupancy relation: x = p/(1−p). This is the Langmuir
relation (Axiom S2e.1), the κ→∞ limit of the master equation. -/
theorem activity_from_occupancy {p x : ℝ} (hp : (1 : ℝ) - p ≠ 0)
    (hx : (1 : ℝ) + x ≠ 0) (h : x / (1 + x) = p) :
    x = p / (1 - p) := by
  field_simp at h ⊢
  linear_combination h

/-- The occupied-site count is binomial at PMF level:
C(Ω,k)·x^k/(1+x)^Ω = C(Ω,k)·p^k·(1−p)^(Ω−k) with p = x/(1+x). -/
theorem occupancy_pmf {Ω k : ℕ} {x p : ℝ} (hx : (1 : ℝ) + x ≠ 0) (hk : k ≤ Ω)
    (hp : p = x / (1 + x)) :
    (Ω.choose k : ℝ) * x ^ k / (1 + x) ^ Ω
      = (Ω.choose k : ℝ) * p ^ k * (1 - p) ^ (Ω - k) := by
  have h1p : (1 : ℝ) - p = 1 / (1 + x) := by
    rw [hp]
    field_simp
  have hpow : ((1 + x)⁻¹) ^ k * ((1 + x)⁻¹) ^ (Ω - k) = ((1 + x)⁻¹) ^ Ω := by
    rw [← pow_add, Nat.add_sub_cancel' hk]
  rw [h1p, hp, div_pow, div_pow, one_pow, one_div, div_eq_mul_inv,
    div_eq_mul_inv, ← inv_pow, ← inv_pow, ← inv_pow]
  conv_lhs => rw [← hpow]
  ring

end TCS
