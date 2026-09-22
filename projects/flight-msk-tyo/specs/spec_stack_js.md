Write ONE JavaScript file `stack.js` (plain ES2015, no modules syntax, no DOM, no external libraries). Output ONLY the code (no markdown fences, no prose). Comments in Russian, short.

## Purpose
Pure helper functions for a stacked spectrum plot.

## Functions (all defined at top level; at the end export them with `if (typeof module !== "undefined") module.exports = { buildStack, russianOrigin };`)

### `buildStack(layers, i0, i1)`
* `layers`: array of `{id, y}` where `y` is an array of numbers (counts per keV per channel, may contain zeros).
* `i0`, `i1`: inclusive channel range used for the contribution.
* `contrib(layer)` = sum of `layer.y[j]` for j from i0 to i1 (skip non-finite values).
* Sort a copy of the layers by contribution in ASCENDING order (smallest first; ties keep input order).
* Return an array of objects `{id, contrib, lower, upper}` in that sorted order, where `upper` is a `Float64Array` of length `layers[0].y.length` holding the running sum of `y` of this layer and of all layers before it in the sorted order, and `lower` is the same array of the previous layer (for the first layer `lower` is a `Float64Array` of zeros of the same length). Non-finite or negative `y` values count as 0.
* Empty input returns `[]`.

### `russianOrigin(label)`
* `label` is a string `"<particle> / <process> / <volume>"` (three parts separated by ` / `). Return a Russian human-readable string `"<particle text>: <process text> · <volume text>"`.
* Particle: `gamma`→`γ-кванты`, `neutron`→`нейтроны`, `e+-`→`электроны и позитроны`, `mu+-`→`мюоны`, `proton`→`протоны`, `other`→`прочие частицы`.
* Process: `primary`→`пришли из воздуха (первичные)`, `eBrem`→`тормозное излучение электронов`, `annihil`→`аннигиляция позитронов`, `compt`→`комптоновское рассеяние`, `conv`→`рождение пар`, `nCapture`→`захват нейтронов`, `RadioactiveDecay`→`радиоактивный распад`, `eIoni`→`ионизация электронами`, `phot`→`фотоэффект`, `protonInelastic`→`неупругое рассеяние протонов`, `neutronInelast`→`неупругое рассеяние нейтронов`, `hadElastic`→`упругое рассеяние адронов`, `other`→`прочие процессы`.
* Volume: `World`→`воздух вне конструкции`, `Pax`→`в слое людей`, `Skin`→`в обшивке`, `Trim`→`в отделке`, `Floor`→`в полу`, `Blanket`→`в изоляции`, `Cargo`→`в багажном отсеке`, `Tube`→`в воздухе у прибора`, `Cabin`→`в воздухе салона`, `other`→`прочее`.
* An unknown word is kept as it is. If the label does not have exactly three parts, return it unchanged.

## Self-test (run when the file is executed by node directly: `if (typeof require !== "undefined" && require.main === module) { ... }`)
* Test 1: layers a=[1,1,1], b=[5,5,5], c=[2,2,2], range 0..2: the order of ids must be a, c, b; `upper` of the last layer must be [8,8,8]; `lower` of the first layer must be [0,0,0]; `upper` of a must be [1,1,1]; `lower` of c must equal `upper` of a.
* Test 2: range 1..1 with layers a=[9,0,0], b=[0,1,0]: the order must be b (contribution 1) after... check that the contributions are `a: 0` and `b: 1`, so the order is a, b.
* Test 3: `russianOrigin("gamma / eBrem / Pax")` must equal `"γ-кванты: тормозное излучение электронов · в слое людей"`; `russianOrigin("x")` must equal `"x"`.
* Print `SELFTEST PASS` and exit with code 0 when all checks hold, else print `SELFTEST FAIL: <which>` and exit with code 1.
