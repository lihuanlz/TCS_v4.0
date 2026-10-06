/-
=====================================================================
S1c: second moment and occupancy variance (finite-size corrections)
=====================================================================

Extends `S1c_PartitionFunction.lean` (first moment) to the second
moment of the occupied-site count N in the grand ensemble of Ω
independent two-state sites. Since each bound ligand occupies exactly
one site, N is simultaneously the bound-ligand number M_L and the
occupied-receptor count Ω_R, so its second moment is the
cross-correlation ⟨M_L Ω_R⟩ behind SI S1c.5, and its variance
Var(N) = Ω·p·(1−p) with p = x/(1+x) is the finite-size fluctuation
behind SI S1c.6:

  - succ_succ_mul_choose: double absorption identity
    (k+2)(k+1)·C(n+2,k+2) = (n+2)(n+1)·C(n,k);
  - secondMoment_eq: Σ_k k²·C(Ω,k)·x^k
      = Ω·x·(1+x)^(Ω−1) + Ω·(Ω−1)·x²·(1+x)^(Ω−2)
    via k² = k(k−1) + k;
  - cross_moment: ⟨M_L Ω_R⟩ = E[N²] = Ω·p + Ω·(Ω−1)·p² at
    p = x/(1+x);
  - occupancy_variance: Var(N) = E[N²] − E[N]² = Ω·p·(1−p),
    the exact binomial variance of the finite system.

Lean 4 + Mathlib, toolchain leanprover/lean4:v4.15.0
-/

import TCS.S1c_PartitionFunction

open Finset

namespace TCS

/-- Second moment of the configuration sum: Σ_k k²·C(Ω,k)·x^k. -/
noncomputable def secondMoment (Ω : ℕ) (x : ℝ) : ℝ :=
  ∑ k ∈ Finset.range (Ω + 1), (k : ℝ) ^ 2 * (Ω.choose k : ℝ) * x ^ k

/-- Double absorption identity:
(k+2)(k+1)·C(n+2,k+2) = (n+2)(n+1)·C(n,k),
two successive applications of Nat.succ_mul_choose_eq. -/
theorem succ_succ_mul_choose (n k : ℕ) :
    ((k + 2 : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ) * (((n + 2).choose (k + 2) : ℕ) : ℝ)
      = ((n + 2 : ℕ) : ℝ) * ((n + 1 : ℕ) : ℝ) * ((n.choose k : ℕ) : ℝ) := by
  have h1 := congrArg (fun m : ℕ => (m : ℝ)) (Nat.succ_mul_choose_eq (n + 1) (k + 1))
  have h2 := congrArg (fun m : ℕ => (m : ℝ)) (Nat.succ_mul_choose_eq n k)
  simp only [Nat.succ_eq_add_one, Nat.cast_mul, Nat.cast_add, Nat.cast_one] at h1 h2
  -- h1 : (k+2)·C(n+2,k+2) = (n+2)·C(n+1,k+1);  h2 : (k+1)·C(n+1,k+1) = (n+1)·C(n,k)
  -- (up to side order); multiply h1 by (k+1) and h2 by (n+2) and add
  linear_combination ((k : ℝ) + 1) * h1 + ((n : ℝ) + 2) * h2

/-- Second-moment sum: Σ_k k²·C(Ω,k)·x^k
= Ω·x·(1+x)^(Ω−1) + Ω·(Ω−1)·x²·(1+x)^(Ω−2).
Split k² = k(k−1) + k; the k part is `occupancyNumerator_eq`, the
k(k−1) part shifts the sum by two and applies the double absorption
identity. -/
theorem secondMoment_eq (Ω : ℕ) (x : ℝ) :
    secondMoment Ω x = (Ω : ℝ) * x * (1 + x) ^ (Ω - 1)
      + (Ω : ℝ) * ((Ω : ℝ) - 1) * x ^ 2 * (1 + x) ^ (Ω - 2) := by
  have hsplit : ∀ k : ℕ, (k : ℝ) ^ 2 * (Ω.choose k : ℝ) * x ^ k
      = ((k : ℝ) * ((k : ℝ) - 1)) * (Ω.choose k : ℝ) * x ^ k
        + (k : ℝ) * (Ω.choose k : ℝ) * x ^ k := by
    intro k
    ring
  rw [secondMoment, Finset.sum_congr rfl (fun k _ => hsplit k), Finset.sum_add_distrib,
    occupancyNumerator_eq Ω x]
  congr 1
  -- Remaining: Σ_k k(k−1)·C(Ω,k)·x^k = Ω·(Ω−1)·x²·(1+x)^(Ω−2)
  cases Ω with
  | zero => simp [Finset.sum_range_succ]
  | succ n =>
    cases n with
    | zero =>
      simp [Finset.sum_range_succ]
    | succ m =>
      -- Ω = m + 2; peel the k = 0 and k = 1 terms (both zero)
      rw [Finset.sum_range_succ', Finset.sum_range_succ']
      simp only [Nat.cast_zero, Nat.cast_one, sub_self, zero_mul, mul_zero,
        zero_sub, add_zero]
      have hterm : ∀ k : ℕ,
          (((k + 2 : ℕ) : ℝ) * (((k + 2 : ℕ) : ℝ) - 1))
            * (((m + 2).choose (k + 2) : ℕ) : ℝ) * x ^ (k + 2)
            = ((m + 2 : ℕ) : ℝ) * ((m + 1 : ℕ) : ℝ) * x ^ 2
              * (((m.choose k : ℕ) : ℝ) * x ^ k) := by
        intro k
        have h := succ_succ_mul_choose m k
        have e1 : ((k + 2 : ℕ) : ℝ) = (k : ℝ) + 2 := by simp
        have e2 : ((k + 2 : ℕ) : ℝ) - 1 = (k : ℝ) + 1 := by simp
        have e3 : ((m + 2 : ℕ) : ℝ) = (m : ℝ) + 2 := by simp
        have e4 : ((m + 1 : ℕ) : ℝ) = (m : ℝ) + 1 := by simp
        rw [e1, e2, e3, e4, pow_add]
        linear_combination (x ^ k * x ^ 2) * h
      rw [Finset.sum_congr rfl (fun k _ => hterm k), ← Finset.mul_sum,
        partitionFunction_factorizes m x]
      have eΩ : m + 1 + 1 - 2 = m := by omega
      have hcast : ((m + 1 + 1 : ℕ) : ℝ) - 1 = ((m + 1 : ℕ) : ℝ) := by simp
      rw [eΩ, hcast]
      ring

/-- Cross-correlation ⟨M_L Ω_R⟩ = E[N²] = Ω·p + Ω·(Ω−1)·p² with
p = x/(1+x) (SI S1c.5): normalized second moment of the occupied-site
count. -/
theorem cross_moment {Ω : ℕ} (hΩ : Ω ≠ 0) {x : ℝ} (hx : (1 : ℝ) + x ≠ 0) :
    secondMoment Ω x / partitionFunction Ω x
      = (Ω : ℝ) * (x / (1 + x))
        + (Ω : ℝ) * ((Ω : ℝ) - 1) * (x / (1 + x)) ^ 2 := by
  rw [secondMoment_eq, partitionFunction_factorizes]
  obtain ⟨m, rfl⟩ := Nat.exists_eq_succ_of_ne_zero hΩ
  have h1 : m + 1 - 1 = m := Nat.add_sub_cancel m 1
  have h2 : ((m + 1 : ℕ) : ℝ) - 1 = (m : ℝ) := by simp
  rw [h1, h2, pow_succ' (1 + x) m]
  cases m with
  | zero =>
    simp only [Nat.cast_zero, pow_zero, mul_one, one_mul, mul_zero, zero_mul,
      add_zero, mul_zero]
    field_simp [hx]
  | succ k =>
    have h3 : k + 1 + 1 - 2 = k := by omega
    rw [h3, pow_succ' (1 + x) k]
    field_simp [hx, pow_ne_zero _ hx]
    ring

/-- Mean occupancy in unnormalized-ratio form:
S₁/Z = Ω·x/(1+x) (equivalent to `mean_occupancy`). -/
theorem mean_occupancy_ratio {Ω : ℕ} (hΩ : Ω ≠ 0) {x : ℝ} (hx : (1 : ℝ) + x ≠ 0) :
    occupancyNumerator Ω x / partitionFunction Ω x = (Ω : ℝ) * (x / (1 + x)) := by
  have hΩr : (Ω : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hΩ
  have hZ : partitionFunction Ω x ≠ 0 := by
    rw [partitionFunction_factorizes]
    exact pow_ne_zero _ hx
  have hm0 := mean_occupancy hΩ hx
  rw [div_eq_iff (mul_ne_zero hΩr hZ)] at hm0
  rw [div_eq_iff hZ, hm0]
  ring

/-- Occupancy variance (SI S1c.6): Var(N) = E[N²] − E[N]²
= Ω·p·(1−p) with p = x/(1+x), the exact finite-size binomial
fluctuation. -/
theorem occupancy_variance {Ω : ℕ} (hΩ : Ω ≠ 0) {x : ℝ} (hx : (1 : ℝ) + x ≠ 0) :
    secondMoment Ω x / partitionFunction Ω x
        - (occupancyNumerator Ω x / partitionFunction Ω x) ^ 2
      = (Ω : ℝ) * (x / (1 + x)) * (1 - x / (1 + x)) := by
  rw [cross_moment hΩ hx, mean_occupancy_ratio hΩ hx]
  ring

end TCS
