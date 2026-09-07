Write two C++17 files for a Geant4 11.2 application: `NpsmLightYield.hh` and
`NpsmLightYield.cc`. Output each as a fenced block preceded by
`=== FILE: <name> ===`. No prose outside the blocks. Comments in Russian,
identifiers in English. Header uses `#pragma once`.

# What this class is

A non-proportional scintillation light-yield model. Given the stopping power S
(MeV/cm) it returns a dimensionless weight in (0, 1]: the fraction of deposited
energy that becomes scintillation light, relative to a proportional detector.

```
                1 - eta * exp[ -(S / S_ons) * exp( -(S_trap / S) ) ]
Weight(S)  =  ------------------------------------------------------
                              1 + S / S_birks
```

CRITICAL, these are the two places implementations get wrong:
- the two exponentials are NESTED: the inner one MULTIPLIES the exponent of the
  outer one. It is NOT `exp(A) * exp(B)`.
- the ratios are INVERTED with respect to each other: in the Onsager term the
  parameter is in the DENOMINATOR (`S / S_ons`); in the trapping term the
  parameter is in the NUMERATOR (`S_trap / S`).

# Class interface

```cpp
class NpsmLightYield {
 public:
  NpsmLightYield();  // parameters default to the NaI(Tl) values below
  // Вес светового выхода. S — тормозная способность, МэВ/см.
  double Weight(double S_MeV_cm) const;
  // Настройка параметров; каждый обязан быть строго положительным,
  // eta — в интервале (0, 1]. Нарушение — громкий отказ, см. ниже.
  void SetParameters(double eta, double sOns, double sTrap, double sBirks);
  // Включена ли модель. Когда выключена, Weight ВСЕГДА возвращает 1.0.
  void SetEnabled(bool on) { fEnabled = on; }
  bool IsEnabled() const { return fEnabled; }
  // Для шапки CSV — постановка обязана лежать В ФАЙЛЕ, а не в его имени.
  double Eta() const; double SOns() const; double STrap() const; double SBirks() const;
 private:
  bool fEnabled = false;   // ВЫКЛЮЧЕНА по умолчанию
  double fEta, fSOns, fSTrap, fSBirks;
};
```

Defaults (maximum a posteriori estimates for NaI(Tl) at 18.8 C):
`eta = 0.596`, `S_ons = 36.4`, `S_trap = 14.6`, `S_birks = 322.0` MeV/cm.

# Behaviour of Weight(S)

1. If `!fEnabled` → return exactly `1.0`. This is the reference branch: with the
   model off the program must reproduce the proportional result bit for bit.
2. Guard the argument. If `S` is not finite, or `S <= 0` → return `1.0` and count
   the occurrence in a mutable counter `mutable long long fNBadS`, exposed by
   `long long BadSCount() const`. Do NOT silently clamp to a small positive
   number: a non-positive stopping power means the caller computed it wrongly,
   and the count must reach the CSV so the defect cannot hide.
   Explain this in a Russian comment.
3. Otherwise evaluate the formula. Compute the inner exponential first, then the
   outer one, exactly in the nested order above.
4. The result is mathematically in (0, 1] for positive parameters. Assert nothing,
   but clamp NOTHING either — if the value leaves that range it means the
   parameters are unphysical and the caller must see it. Instead count values
   outside `(0.0, 1.0]` in a second mutable counter `fNOutOfRange`, exposed by
   `long long OutOfRangeCount() const`, and return the value as computed.

# SetParameters validation

Reject with a loud failure — `std::fprintf(stderr, "NpsmLightYield: FATAL ...")`
followed by `std::abort()` — when any of:
- `eta <= 0.0` or `eta > 1.0`;
- `sOns <= 0.0`, `sTrap <= 0.0`, `sBirks <= 0.0`;
- any argument is not finite (`std::isfinite`).
The message must name the offending parameter and its value. A Russian comment
explains why abort rather than a return code: a physics run started with
nonsensical parameters produces a plausible-looking spectrum, and a silently
wrong spectrum is worse than no spectrum.

# Numerical notes to implement

- Use `std::exp`, `std::isfinite` from `<cmath>`.
- For very large `S` the outer exponent tends to 0 and the numerator tends to
  `1 - eta`; for very small `S` the inner exponential underflows to 0, the outer
  exponent becomes 0 and the numerator also tends to `1 - eta`. Both limits are
  finite — no special-casing needed, but state this in a comment so a reader does
  not add a spurious guard.
- The denominator `1 + S/S_birks` is always >= 1 for positive S, so no division
  by zero is possible once S > 0 is established.

# Self-test entry point

Add to the .cc a function with C linkage-free signature:

```cpp
// Автопроверка формы кривой. Возвращает 0, если все утверждения верны,
// иначе номер первого нарушенного. Печатает подробности в stdout.
int NpsmLightYieldSelfTest();
```

It must check, with the default NaI(Tl) parameters and the model ENABLED:
1. `Weight(1.0)` and `Weight(1000.0)` are both finite and in `(0, 1]`;
2. `Weight` is strictly decreasing between `S = 100` and `S = 1000`
   (higher stopping power → stronger Birks quenching → less light);
3. with `eta` set to a value near zero (use `SetParameters(1e-9, 36.4, 14.6, 322.0)`)
   the numerator is essentially 1, so `Weight(S)` approaches `1/(1 + S/S_birks)`
   within 1e-6 at `S = 200`;
4. with the model DISABLED, `Weight(12345.0)` equals exactly `1.0`;
5. `Weight(-1.0)` returns `1.0` and increments the bad-S counter by exactly one.
Print each check as `selftest N: OK` or `selftest N: FAIL <details>`.
