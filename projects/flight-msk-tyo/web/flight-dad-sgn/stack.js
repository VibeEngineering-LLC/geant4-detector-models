function buildStack(layers, i0, i1) {
  if (layers.length === 0) return [];

  const copy = layers.map(layer => ({
    id: layer.id,
    y: layer.y,
    contrib: 0
  }));

  for (let i = 0; i < copy.length; i++) {
    let sum = 0;
    for (let j = i0; j <= i1; j++) {
      const val = copy[i].y[j];
      if (Number.isFinite(val) && val >= 0) {
        sum += val;
      }
    }
    copy[i].contrib = sum;
  }

  copy.sort((a, b) => {
    if (a.contrib === b.contrib) return 0;
    return a.contrib < b.contrib ? -1 : 1;
  });

  const result = [];
  let runningSum = new Float64Array(copy[0].y.length).fill(0);

  for (let i = 0; i < copy.length; i++) {
    const lower = i === 0 ? new Float64Array(runningSum.length).fill(0) : runningSum.slice();
    const upper = new Float64Array(runningSum.length);
    
    for (let j = 0; j < runningSum.length; j++) {
      const val = copy[i].y[j];
      const v = Number.isFinite(val) && val >= 0 ? val : 0;
      runningSum[j] += v;
      upper[j] = runningSum[j];
    }

    result.push({
      id: copy[i].id,
      contrib: copy[i].contrib,
      lower: lower,
      upper: upper
    });
  }

  return result;
}

const ORIGIN_DICT_RU = {
  particle: { gamma: "γ-кванты", neutron: "нейтроны", "e+-": "электроны и позитроны", "mu+-": "мюоны", proton: "протоны", other: "прочие частицы" },
  process: { primary: "пришли из воздуха (первичные)", eBrem: "тормозное излучение электронов", annihil: "аннигиляция позитронов",
    compt: "комптоновское рассеяние", conv: "рождение пар", nCapture: "захват нейтронов", RadioactiveDecay: "радиоактивный распад",
    eIoni: "ионизация электронами", phot: "фотоэффект", protonInelastic: "неупругое рассеяние протонов",
    neutronInelastic: "неупругое рассеяние нейтронов", hadElastic: "упругое рассеяние адронов",
    muPairProd: "рождение пар мюоном", photonNuclear: "фотоядерная реакция", muIoni: "ионизация мюонами",
    muBrems: "тормозное излучение мюонов", other: "прочие процессы" },
  volume: { World: "воздух вне конструкции", Pax: "в слое людей", Skin: "в обшивке", Trim: "в отделке", Floor: "в полу",
    Blanket: "в изоляции", Cargo: "в багажном отсеке", Tube: "в воздухе у прибора", Cabin: "в воздухе салона", other: "прочее" }
};

const ORIGIN_DICT_EN = {
  particle: { gamma: "γ-rays", neutron: "neutrons", "e+-": "electrons and positrons", "mu+-": "muons", proton: "protons", other: "other particles" },
  process: { primary: "arrived from outside air (primary)", eBrem: "electron bremsstrahlung", annihil: "positron annihilation",
    compt: "Compton scattering", conv: "pair production", nCapture: "neutron capture", RadioactiveDecay: "radioactive decay",
    eIoni: "electron ionization", phot: "photoelectric effect", protonInelastic: "proton inelastic scattering",
    neutronInelastic: "neutron inelastic scattering", hadElastic: "hadron elastic scattering",
    muPairProd: "muon pair production", photonNuclear: "photonuclear reaction", muIoni: "muon ionization",
    muBrems: "muon bremsstrahlung", other: "other processes" },
  volume: { World: "air outside the airframe", Pax: "in the passenger layer", Skin: "in the skin", Trim: "in the trim", Floor: "in the floor",
    Blanket: "in the insulation", Cargo: "in the cargo hold", Tube: "in the air near the instrument", Cabin: "in the cabin air", other: "other" }
};

function russianOrigin(label, lang) {
  const parts = label.split(" / ");
  if (parts.length !== 3) return label;
  const [particle, process, volume] = parts;
  const dict = lang === "en" ? ORIGIN_DICT_EN : ORIGIN_DICT_RU;
  const p = dict.particle[particle] || particle;
  const pr = dict.process[process] || process;
  const v = dict.volume[volume] || volume;
  return `${p}: ${pr} · ${v}`;
}

if (typeof require !== "undefined" && require.main === module) {
  // Test 1
  const layers1 = [
    { id: "a", y: [1, 1, 1] },
    { id: "b", y: [5, 5, 5] },
    { id: "c", y: [2, 2, 2] }
  ];
  const res1 = buildStack(layers1, 0, 2);
  if (res1[0].id !== "a" || res1[1].id !== "c" || res1[2].id !== "b") {
    console.log("SELFTEST FAIL: Test 1 order");
    process.exit(1);
  }
  if (res1[2].upper[0] !== 8) {
    console.log("SELFTEST FAIL: Test 1 upper");
    process.exit(1);
  }
  if (res1[0].lower[0] !== 0) {
    console.log("SELFTEST FAIL: Test 1 lower first");
    process.exit(1);
  }
  if (res1[0].upper[0] !== 1) {
    console.log("SELFTEST FAIL: Test 1 upper a");
    process.exit(1);
  }
  if (!res1[1].lower.every((v, i) => v === res1[0].upper[i])) {
    console.log("SELFTEST FAIL: Test 1 lower c");
    process.exit(1);
  }

  // Test 2
  const layers2 = [
    { id: "a", y: [9, 0, 0] },
    { id: "b", y: [0, 1, 0] }
  ];
  const res2 = buildStack(layers2, 1, 1);
  if (res2[0].id !== "a" || res2[1].id !== "b") {
    console.log("SELFTEST FAIL: Test 2 order");
    process.exit(1);
  }

  // Test 3
  const r1 = russianOrigin("gamma / eBrem / Pax");
  if (r1 !== "γ-кванты: тормозное излучение электронов · в слое людей") {
    console.log("SELFTEST FAIL: Test 3 russianOrigin");
    process.exit(1);
  }
  const r2 = russianOrigin("x");
  if (r2 !== "x") {
    console.log("SELFTEST FAIL: Test 3 russianOrigin unchanged");
    process.exit(1);
  }
  const r3 = russianOrigin("gamma / eBrem / Pax", "en");
  if (r3 !== "γ-rays: electron bremsstrahlung · in the passenger layer") {
    console.log("SELFTEST FAIL: Test 4 russianOrigin en");
    process.exit(1);
  }

  console.log("SELFTEST PASS");
  process.exit(0);
}

if (typeof module !== "undefined") module.exports = { buildStack, russianOrigin };
