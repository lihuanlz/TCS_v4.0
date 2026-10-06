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
  have h1 : ((n + 2 : ℕ) : ℝ) * (((n + 1).choose (k + 1) : ℕ) : ℝ)
      = (((n + 2).choose (k + 2) : ℕ) : ℝ) * ((k + 2 : ℕ) : ℝ) := by
    have h := congrArg (fun m : ℕ => (m : ℝ)) (Nat.succ_mul_choose_eq (n + 1) (k + 1))
    simp only [Nat.succ_eq_add_one, Nat.cast_mul] at h
    exact h
  have h2 : ((n + 1 : ℕ) : ℝ) * ((n.choose k : ℕ) : ℝ)
      = (((n + 1).choose (k + 1) : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ) := by
    have h := congrArg (fun m : ℕ => (m : ℝ)) (Nat.succ_mul_choose_eq n k)
    simp only [Nat.succ_eq_add_one, Nat.cast_mul] at h
    exact h
  calc ((k + 2 : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ) * (((n + 2).choose (k + 2) : ℕ) : ℝ)
      = ((k + 1 : ℕ) : ℝ)
          * ((((n + 2).choose (k + 2) : ℕ) : ℝ) * ((k + 2 : ℕ) : ℝ)) := by ring
    _ = ((k + 1 : ℕ) : ℝ)
          * (((n + 2 : ℕ) : ℝ) * (((n + 1).choose (k + 1) : ℕ) : ℝ)) := by rw [← h1]
    _ = ((n + 2 : ℕ) : ℝ)
          * ((((n + 1).choose (k + 1) : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ)) := by ring
    _ = ((n + 2 : ℕ) : ℝ) * (((n + 1 : ℕ) : ℝ) * ((n.choose k : ℕ) : ℝ)) := by rw [← h2]
    _ = ((n + 2 : ℕ) : ℝ) * ((n + 1 : ℕ) : ℝ) * ((n.choose k : ℕ) : ℝ) := by ring

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
  rw [secondMoment, Finset.sum_congr rfl (fun k _ => hsplit k), Finset.sum_add_distrib]
  rw [show (∑ k ∈ Finset.range (Ω + 1), (k : ℝ) * (Ω.choose k : ℝ) * x ^ k)
        = occupancyNumerator Ω x from rfl, occupancyNumerator_eq Ω x]
  -- Goal: S + T = T + D with S the k(k−1) sum, T = Ω·x·(1+x)^(Ω−1),
  -- D = Ω·(Ω−1)·x²·(1+x)^(Ω−2).  We prove S = D per case and close
  -- by commutativity.
  cases Ω with
  | zero => simp [Finset.sum_range_succ]
  | succ n =>
    cases n with
    | zero => simp [Finset.sum_range_succ]
    | succ m =>
      -- Normalize the goal to (m+2)-form so all cast atoms match.
      show (∑ k ∈ Finset.range (m + 2 + 1), (k : ℝ) * ((k : ℝ) - 1)
              * ((m + 2).choose k : ℝ) * x ^ k)
            + ((m + 2 : ℕ) : ℝ) * x * (1 + x) ^ (m + 2 - 1)
          = ((m + 2 : ℕ) : ℝ) * x * (1 + x) ^ (m + 2 - 1)
            + ((m + 2 : ℕ) : ℝ) * (((m + 2 : ℕ) : ℝ) - 1) * x ^ 2 * (1 + x) ^ (m + 2 - 2)
      -- Peel the k = 0 and k = 1 terms (both vanish since k(k−1) = 0).
      rw [Finset.sum_range_succ', Finset.sum_range_succ']
      have hf0 : ((0 : ℕ) : ℝ) * (((0 : ℕ) : ℝ) - 1)
          * (((m + 2).choose 0 : ℕ) : ℝ) * x ^ 0 = 0 := by simp
      have hf1 : ((0 + 1 : ℕ) : ℝ) * (((0 + 1 : ℕ) : ℝ) - 1)
          * (((m + 2).choose (0 + 1) : ℕ) : ℝ) * x ^ (0 + 1) = 0 := by simp
      rw [hf0, hf1]
      simp only [add_zero]
      have hterm : ∀ k : ℕ,
          (((k + 2 : ℕ) : ℝ) * (((k + 2 : ℕ) : ℝ) - 1))
            * (((m + 2).choose (k + 2) : ℕ) : ℝ) * x ^ (k + 2)
            = ((m + 2 : ℕ) : ℝ) * ((m + 1 : ℕ) : ℝ) * x ^ 2
              * (((m.choose k : ℕ) : ℝ) * x ^ k) := by
        intro k
        have h := succ_succ_mul_choose m k
        have e2 : ((k + 2 : ℕ) : ℝ) - 1 = ((k + 1 : ℕ) : ℝ) := by
          simp only [Nat.cast_add, Nat.cast_one, Nat.cast_ofNat]; ring
        rw [e2, pow_add]
        linear_combination (x ^ k * x ^ 2) * h
      rw [Finset.sum_congr rfl (fun k _ => hterm k), ← Finset.mul_sum,
        show (∑ i ∈ Finset.range (m + 1), ((m.choose i : ℕ) : ℝ) * x ^ i)
          = partitionFunction m x from rfl, partitionFunction_factorizes m x]
      have eΩ : m + 2 - 2 = m := by omega
      have hcast : ((m + 2 : ℕ) : ℝ) - 1 = ((m + 1 : ℕ) : ℝ) := by
        simp only [Nat.cast_add, Nat.cast_one, Nat.cast_ofNat]; ring
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
