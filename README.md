# Vapor–Liquid Equilibrium: Methanol–Water

Calculates bubble point temperatures and equilibrium vapor compositions for a
binary methanol–water mixture at constant pressure, using Raoult's law with
Antoine vapor pressure correlations.

## The problem

For a binary mixture at fixed total pressure `P`, the bubble point temperature
is the temperature at which the partial pressures of both components first sum
to `P`. Since the vapor mole fractions must add to one, this gives the bubble
point condition:

```
x₁·P₁°(T) + (1 - x₁)·P₂°(T) - P = 0
```

This equation is implicit in `T` and has no closed-form solution: the unknown.

It is solved numerically with Brent's method. The search interval is bounded by
the normal boiling points of the two pure components — the mixture must boil
somewhere between them — which means the root is guaranteed to be bracketed and
no initial guess is required.

Once `T` is known, the vapor composition follows directly from Raoult's law:

```
y₁ = x₁·P₁°(T) / P
```

## Results

At 760 mmHg:

| x₁ (liquid) | T (°C) | y₁ (vapor) |
|------------:|-------:|-----------:|
| 0.00 | 100.00 | 0.0000 |
| 0.10 |  93.70 | 0.2849 |
| 0.30 |  84.02 | 0.6164 |
| 0.50 |  76.82 | 0.7953 |
| 0.70 |  71.17 | 0.9032 |
| 0.90 |  66.56 | 0.9736 |
| 1.00 |  64.55 | 1.0000 |

A liquid containing 50 mol% methanol boils at 76.8 °C and is in equilibrium with
a vapor containing 79.5 mol% methanol. That enrichment is what makes distillation
work, and the vertical distance between the equilibrium curve and the `y = x`
diagonal is the separation achieved by one theoretical stage.

The gap is widest at low methanol concentrations and narrows sharply near
purity: going from x₁ = 0.10 to y₁ = 0.28 in a single stage is nearly a
tripling, while x₁ = 0.90 yields only y₁ = 0.97. This is why distillation
columns need progressively more stages toward the top.

![x-y diagram](Equilibrium_Curve.png)

![T-x-y diagram](Bubble_Point_Temperature.png)

## Assumptions and limitations

The calculation assumes an ideal liquid phase (Raoult's law) and an ideal vapor
phase. Methanol–water is not an ideal system — it shows positive deviation from
Raoult's law due to hydrogen bonding — so the predicted vapor is less enriched
in methanol than experiment shows, and the predicted bubble point temperatures
are slightly high.

## Usage

```bash
pip install -r requirements.txt
python vle.py
```

Prints a bubble point table and writes both diagrams as PNG files.

## Data

Antoine coefficients, valid for `P` in mmHg and `T` in °C, in the form
`log₁₀(P°) = A - B/(T + C)`:

| Component | A | B | C |
|-----------|---------|----------|---------|
| Methanol  | 8.08097 | 1582.271 | 239.726 |
| Water     | 8.07131 | 1730.630 | 233.426 |

*(Add the source you took these from — a textbook appendix, NIST, or Perry's
Handbook. Antoine coefficients vary between sources and are only valid over a
limited temperature range, so the citation matters.)*

## Verification

The code checks itself on three things before computing anything:

- `antoine_inverse(100, water)` returns 760 mmHg — the known normal boiling point
- `antoine(760, methanol)` returns 64.55 °C — literature value 64.7 °C
- `antoine(antoine_inverse(T))` returns `T` — the two Antoine forms are mutual inverses

At the endpoints, `x₁ = 0` and `x₁ = 1` return `y₁ = 0` and `y₁ = 1` exactly,
as they must for pure components.
