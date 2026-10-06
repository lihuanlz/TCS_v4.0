/-
=====================================================================
S2e.10.2: the Jacobian determinant condition (Eq. S2e.10b)
=====================================================================

This module formalizes the determinant condition that
`S2e10_Temperature.lean` previously declared "deliberately NOT
formalized": for the digital counting readout at three temperatures,
the log-parameter Jacobian of the observables

  Y_j = M / (1 + κ₀·g_j),   g_j = exp(h·x_j),
  x_j = −(1/R)(1/T_j − 1/T₀),  h = ΔH,

with respect to (ln M, ln κ₀, h) has determinant

  det J = Y₁·Y₂·Y₃ · D,   D = u₁u₂(x₂−x₁) + u₂u₃(x₃−x₂) + u₁u₃(x₁−x₃),

where u_j = κ₀g_j/(1+κ₀g_j) is the depletion share. The scientific
content of Eq. (S2e.10b) is:

  * `jacobian_expansion`: the 3×3 cofactor expansion equals
    Y₁Y₂Y₃·D (pure algebra);
  * `tempJacDet_zero_of_dH_zero`: ΔH = 0 ⇒ all u_j coincide ⇒ D = 0
    (the degeneracy persists);
  * `tempJacDet_ne_zero`: κ₀ > 0, h ≠ 0 and pairwise distinct
    temperatures ⟹ D ≠ 0;
  * `jacobian_full_ne_zero`: the full determinant is nonzero whenever
    M ≠ 0 as well.

The nonvanishing proof is Rolle's theorem applied twice to
N(x) = α + κ₀(α+β+γx)e^{hx}: if D = 0, the explicit cofactor triple
(α, β, γ) = (u₂u₃(x₃−x₂), u₂x₂−u₃x₃, u₃−u₂) makes N vanish at all three
temperatures; N'(x) = κ₀e^{hx}·(γ + h(α+β) + hγx) is an exponential
times a *linear* function, which cannot vanish at two distinct Rolle
points unless γ = 0; but γ = u₃ − u₂ ≠ 0 because u is injective.
Hence three temperatures with distinct reduced reciprocals make
(M, κ₀, ΔH) structurally identifiable: the SI claim "nonzero whenever
ΔH ≠ 0 and the temperatures are distinct".

Lean 4 + Mathlib, toolchain leanprover/lean4:v4.15.0
-/

import TCS.S2e10_Temperature

namespace TCS

/-- Depletion share at reduced reciprocal temperature x:
u = κ₀·e^(hx)/(1 + κ₀·e^(hx)), the fraction of captured targets at a
single temperature. This is the S2e.10.2 analogue of the digital
depletion factor of Theorem S2c.4.1. -/
noncomputable def deplU (κ₀ h x : ℝ) : ℝ :=
  κ₀ * Real.exp (h * x) / (1 + κ₀ * Real.exp (h * x))

/-- The reduced log-parameter Jacobian determinant (the full Jacobian
carries the additional factor Y₁Y₂Y₃, see `jacobian_expansion`):
D = u₁u₂(x₂−x₁) + u₂u₃(x₃−x₂) + u₁u₃(x₁−x₃). -/
def tempJacDet (u₁ u₂ u₃ x₁ x₂ x₃ : ℝ) : ℝ :=
  u₁ * u₂ * (x₂ - x₁) + u₂ * u₃ * (x₃ - x₂) + u₁ * u₃ * (x₁ - x₃)

/-- The depletion share is injective in the temperature variable when
κ₀ > 0 and ΔH ≠ 0: equal shares force equal reduced temperatures,
because u = v/(1+v) is injective in v and v = κ₀e^{hx} is injective
in x. Pure algebra plus injectivity of exp. -/
theorem deplU_injective {κ₀ h : ℝ} (hκ : 0 < κ₀) (hh : h ≠ 0) :
    Function.Injective (deplU κ₀ h) := by
  intro a b hab
  simp only [deplU] at hab
  have hpa : (0 : ℝ) < 1 + κ₀ * Real.exp (h * a) := by positivity
  have hpb : (0 : ℝ) < 1 + κ₀ * Real.exp (h * b) := by positivity
  field_simp [ne_of_gt hpa, ne_of_gt hpb] at hab
  have h2 : κ₀ * Real.exp (h * a) = κ₀ * Real.exp (h * b) := by
    linear_combination hab
  have h3 : Real.exp (h * a) = Real.exp (h * b) :=
    mul_left_cancel₀ (ne_of_gt hκ) h2
  have h4 : h * a = h * b := Real.exp_injective h3
  exact mul_left_cancel₀ hh h4

/-- S2e.10.2 degeneracy direction: at ΔH = 0 every temperature factor
is 1, all depletion shares coincide, and the determinant vanishes —
no temperature protocol helps when the affinity is
temperature-independent. -/
theorem tempJacDet_zero_of_dH_zero (κ₀ : ℝ) (x₁ x₂ x₃ : ℝ) :
    tempJacDet (deplU κ₀ 0 x₁) (deplU κ₀ 0 x₂) (deplU κ₀ 0 x₃) x₁ x₂ x₃ = 0 := by
  simp only [deplU, zero_mul, Real.exp_zero, mul_one, tempJacDet]
  ring

/-- Determinant expansion: the full log-parameter Jacobian determinant
equals Y₁Y₂Y₃ times the reduced determinant D. The rows of J are
(Y_j, −Y_j·u_j, −Y_j·u_j·x_j); cofactor expansion along the first
column and pulling out the row factors gives this identity
(pure algebra). -/
theorem jacobian_expansion (Y₁ Y₂ Y₃ u₁ u₂ u₃ x₁ x₂ x₃ : ℝ) :
    Y₁ * Y₂ * Y₃ * tempJacDet u₁ u₂ u₃ x₁ x₂ x₃
      = Y₁ * ((-Y₂ * u₂) * (-Y₃ * u₃ * x₃) - (-Y₃ * u₃) * (-Y₂ * u₂ * x₂))
        - Y₂ * ((-Y₁ * u₁) * (-Y₃ * u₃ * x₃) - (-Y₃ * u₃) * (-Y₁ * u₁ * x₁))
        + Y₃ * ((-Y₁ * u₁) * (-Y₂ * u₂ * x₂) - (-Y₂ * u₂) * (-Y₁ * u₁ * x₁)) := by
  simp only [tempJacDet]
  ring

/-- Rolle core: if N(x) = α + κ₀(α+β+γx)e^{hx} vanishes at three
strictly increasing points and κ₀, h ≠ 0, then γ = 0. Between each
pair of zeros, Rolle's theorem gives a zero of
N'(x) = κ₀e^{hx}·(γ + h(α+β) + hγx); the bracket is linear in x with
slope hγ, and two distinct zeros of a linear function force the slope
to vanish. -/
theorem rolle_linear_vanishing {α β γ κ₀ h : ℝ} (hκ : κ₀ ≠ 0) (hh : h ≠ 0)
    {x₁ x₂ x₃ : ℝ} (hx12 : x₁ < x₂) (hx23 : x₂ < x₃)
    (h1 : α + κ₀ * ((α + β + γ * x₁) * Real.exp (h * x₁)) = 0)
    (h2 : α + κ₀ * ((α + β + γ * x₂) * Real.exp (h * x₂)) = 0)
    (h3 : α + κ₀ * ((α + β + γ * x₃) * Real.exp (h * x₃)) = 0) :
    γ = 0 := by
  set N : ℝ → ℝ := fun y => α + κ₀ * ((α + β + γ * y) * Real.exp (h * y)) with hN
  have hderiv : ∀ x : ℝ,
      deriv N x = κ₀ * Real.exp (h * x) * (γ + h * (α + β + γ * x)) := by
    intro x
    have ha : HasDerivAt (fun y : ℝ => α + β + γ * y) (γ * 1) x :=
      ((hasDerivAt_id x).const_mul γ).const_add (α + β)
    have hb : HasDerivAt (fun y : ℝ => Real.exp (h * y))
        (Real.exp (h * x) * (h * 1)) x :=
      (Real.hasDerivAt_exp (h * x)).comp x ((hasDerivAt_id x).const_mul h)
    have hc := (ha.mul hb).const_mul κ₀
    have hd := hc.const_add α
    rw [hN, hd.deriv]
    ring
  obtain ⟨ξ₁, ⟨hξ₁lo, hξ₁hi⟩, hξ₁0⟩ :=
    exists_deriv_eq_zero (f := N) hx12 (by rw [hN]; fun_prop)
      (by rw [show N x₁ = 0 from h1, show N x₂ = 0 from h2])
  obtain ⟨ξ₂, ⟨hξ₂lo, hξ₂hi⟩, hξ₂0⟩ :=
    exists_deriv_eq_zero (f := N) hx23 (by rw [hN]; fun_prop)
      (by rw [show N x₂ = 0 from h2, show N x₃ = 0 from h3])
  rw [hderiv ξ₁] at hξ₁0
  rw [hderiv ξ₂] at hξ₂0
  have hL1 : γ + h * (α + β + γ * ξ₁) = 0 := by
    rcases mul_eq_zero.mp hξ₁0 with hκe | hL
    · exact absurd hκe (mul_ne_zero hκ (Real.exp_ne_zero _))
    · exact hL
  have hL2 : γ + h * (α + β + γ * ξ₂) = 0 := by
    rcases mul_eq_zero.mp hξ₂0 with hκe | hL
    · exact absurd hκe (mul_ne_zero hκ (Real.exp_ne_zero _))
    · exact hL
  have hmid : h * (α + β + γ * ξ₁) = h * (α + β + γ * ξ₂) := by linarith
  have hkey : h * γ * (ξ₂ - ξ₁) = 0 := by linear_combination -hmid
  have hξ : ξ₂ - ξ₁ ≠ 0 := sub_ne_zero.mpr (ne_of_gt (lt_trans hξ₁hi hξ₂lo))
  rcases mul_eq_zero.mp hkey with hhγ | hξ'
  · exact (mul_eq_zero.mp hhγ).resolve_left hh
  · exact absurd hξ' hξ

/-- S2e.10.2 (Eq. S2e.10b), main result: for κ₀ > 0, ΔH ≠ 0 and three
strictly increasing (reduced reciprocal) temperatures, the reduced
Jacobian determinant is nonzero — three temperatures make
(M, κ₀, ΔH) structurally identifiable from digital counting readout. -/
theorem tempJacDet_ne_zero {κ₀ h : ℝ} (hκ : 0 < κ₀) (hh : h ≠ 0)
    {x₁ x₂ x₃ : ℝ} (hx12 : x₁ < x₂) (hx23 : x₂ < x₃) :
    tempJacDet (deplU κ₀ h x₁) (deplU κ₀ h x₂) (deplU κ₀ h x₃) x₁ x₂ x₃ ≠ 0 := by
  intro HD
  have hγ0 : deplU κ₀ h x₃ - deplU κ₀ h x₂ ≠ 0 :=
    sub_ne_zero.mpr ((deplU_injective hκ hh).ne (ne_of_gt hx23))
  have hF : ∀ y : ℝ, ∀ α β γ : ℝ,
      (α + β * deplU κ₀ h y + γ * deplU κ₀ h y * y) * (1 + κ₀ * Real.exp (h * y))
        = α + κ₀ * ((α + β + γ * y) * Real.exp (h * y)) := by
    intro y α β γ
    have hpos : (0 : ℝ) < 1 + κ₀ * Real.exp (h * y) := by positivity
    simp only [deplU]
    field_simp [ne_of_gt hpos]
    ring
  set α : ℝ := deplU κ₀ h x₂ * deplU κ₀ h x₃ * (x₃ - x₂) with hα
  set β : ℝ := deplU κ₀ h x₂ * x₂ - deplU κ₀ h x₃ * x₃ with hβ
  set γ : ℝ := deplU κ₀ h x₃ - deplU κ₀ h x₂ with hγ
  have hN1 : α + κ₀ * ((α + β + γ * x₁) * Real.exp (h * x₁)) = 0 := by
    rw [← hF x₁ α β γ]
    have hid : α + β * deplU κ₀ h x₁ + γ * deplU κ₀ h x₁ * x₁
        = tempJacDet (deplU κ₀ h x₁) (deplU κ₀ h x₂) (deplU κ₀ h x₃) x₁ x₂ x₃ := by
      rw [hα, hβ, hγ]
      unfold tempJacDet
      ring
    rw [hid, HD]
    ring
  have hN2 : α + κ₀ * ((α + β + γ * x₂) * Real.exp (h * x₂)) = 0 := by
    rw [← hF x₂ α β γ]
    have hid : α + β * deplU κ₀ h x₂ + γ * deplU κ₀ h x₂ * x₂ = 0 := by
      rw [hα, hβ, hγ]
      ring
    rw [hid]
    ring
  have hN3 : α + κ₀ * ((α + β + γ * x₃) * Real.exp (h * x₃)) = 0 := by
    rw [← hF x₃ α β γ]
    have hid : α + β * deplU κ₀ h x₃ + γ * deplU κ₀ h x₃ * x₃ = 0 := by
      rw [hα, hβ, hγ]
      ring
    rw [hid]
    ring
  exact hγ0 (rolle_linear_vanishing (ne_of_gt hκ) hh hx12 hx23 hN1 hN2 hN3)

/-- S2e.10.2, full form: the log-parameter Jacobian determinant of the
three-temperature digital observables Y_j = M/(1+κ₀g_j) is nonzero
whenever M ≠ 0, κ₀ > 0, ΔH ≠ 0 and the temperatures are distinct.
This is exactly the nonzero claim of Eq. (S2e.10b). -/
theorem jacobian_full_ne_zero {M κ₀ h : ℝ} (hM : M ≠ 0) (hκ : 0 < κ₀) (hh : h ≠ 0)
    {x₁ x₂ x₃ : ℝ} (hx12 : x₁ < x₂) (hx23 : x₂ < x₃) :
    Ydig M κ₀ (Real.exp (h * x₁)) * Ydig M κ₀ (Real.exp (h * x₂))
        * Ydig M κ₀ (Real.exp (h * x₃))
        * tempJacDet (deplU κ₀ h x₁) (deplU κ₀ h x₂) (deplU κ₀ h x₃) x₁ x₂ x₃
      ≠ 0 := by
  have hY : ∀ x : ℝ, Ydig M κ₀ (Real.exp (h * x)) ≠ 0 := by
    intro x
    have hpos : (0 : ℝ) < 1 + κ₀ * Real.exp (h * x) := by positivity
    simp only [Ydig]
    exact div_ne_zero hM (ne_of_gt hpos)
  exact mul_ne_zero
    (mul_ne_zero (mul_ne_zero (hY x₁) (hY x₂)) (hY x₃))
    (tempJacDet_ne_zero hκ hh hx12 hx23)

end TCS
