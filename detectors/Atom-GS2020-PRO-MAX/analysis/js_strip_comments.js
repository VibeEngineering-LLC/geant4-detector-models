// #GS-79: снимает комментарии из JS парсером TypeScript (removeComments), код и строки не трогает.
// Запуск: node js_strip_comments.js <in.js> <out.js>; TS — глобальный npm (typescript@6).
const fs = require("fs"), path = require("path");
const ts = require(path.join(process.env.APPDATA, "npm", "node_modules", "typescript"));
const src = fs.readFileSync(process.argv[2], "utf8");
const r = ts.transpileModule(src, { compilerOptions: { removeComments: true, target: ts.ScriptTarget.ESNext,
  module: ts.ModuleKind.ESNext, allowJs: true, alwaysStrict: false, strict: false, ignoreDeprecations: "6.0", noEmitHelpers: true, newLine: ts.NewLineKind.LineFeed }, reportDiagnostics: true });
const errs = (r.diagnostics || []).filter(d => d.category === ts.DiagnosticCategory.Error);
if (errs.length) { console.error("ОТКАЗ: синтаксис " + errs.map(d => ts.flattenDiagnosticMessageText(d.messageText, " ")).join("; ")); process.exit(2); }
fs.writeFileSync(process.argv[3], r.outputText.replace(/^export \{\};\s*$/m, ""), "utf8");
