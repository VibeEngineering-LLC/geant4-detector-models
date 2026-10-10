// #GS-79: проверка снятия комментариев — синтаксические деревья двух JS совпадают (вид узла + текст листьев).
// Запуск: node js_ast_equal.js <a.js> <b.js>; печатает число узлов и первое расхождение; код 1 при расхождении.
const fs = require("fs"), path = require("path");
const ts = require(path.join(process.env.APPDATA, "npm", "node_modules", "typescript"));
function flat(fn) {
  const src = fs.readFileSync(fn, "utf8"), sf = ts.createSourceFile(fn, src, ts.ScriptTarget.ESNext, true, ts.ScriptKind.JS), out = [];
  (function walk(n) {
    const kids = n.getChildren(sf);
    out.push(kids.length ? ts.SyntaxKind[n.kind] : ts.SyntaxKind[n.kind] + ":" + src.slice(n.getStart(sf), n.end));
    kids.forEach(walk);
  })(sf);
  return out;
}
const a = flat(process.argv[2]), b = flat(process.argv[3]);
let i = 0; while (i < a.length && i < b.length && a[i] === b[i]) i++;
if (i === a.length && i === b.length) { console.log("равны: узлов " + a.length); process.exit(0); }
console.log("РАСХОЖДЕНИЕ на узле " + i + " из " + a.length + "/" + b.length + ": " + a[i] + " | " + b[i]); process.exit(1);
