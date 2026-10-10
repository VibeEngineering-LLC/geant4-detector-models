// #GS-79: литералы JS (строки, шаблоны, regex) с кириллицей — сырой текст литерала и число вхождений по файлу.
// Запуск: node js_cyr_literals.js <out.json> <a.js> [b.js ...]; парсер TypeScript (глобальный npm).
const fs = require("fs"), path = require("path");
const ts = require(path.join(process.env.APPDATA, "npm", "node_modules", "typescript"));
const CYR = /[А-Яа-яЁё]/, K = ts.SyntaxKind, out = {};
for (const fn of process.argv.slice(3)) {
  const src = fs.readFileSync(fn, "utf8"), sf = ts.createSourceFile(fn, src, ts.ScriptTarget.ESNext, true, ts.ScriptKind.JS);
  const lits = {};
  (function walk(n) {
    if ([K.StringLiteral, K.NoSubstitutionTemplateLiteral, K.TemplateHead, K.TemplateMiddle, K.TemplateTail, K.RegularExpressionLiteral].includes(n.kind)) {
      const raw = src.slice(n.getStart(sf), n.end);
      if (CYR.test(raw)) lits[raw] = (lits[raw] || 0) + 1;
    }
    ts.forEachChild(n, walk);
  })(sf);
  out[path.basename(fn)] = lits;
}
fs.writeFileSync(process.argv[2], JSON.stringify(out, null, 1), "utf8");
