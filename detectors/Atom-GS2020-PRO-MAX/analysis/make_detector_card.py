# #DET-1: генератор карточки детектора GS2020 PRO MAX + Маринелли 1 л (РАДЭК) из референса тория (ОИСН-10, 26.09.2026).
# Значения — таблица «Параметры пиков» СпектраЛайн (README-референсы.md, скрин 2026-09-26) × множители #SHAPE-1 (gs2020_fit_env.sh).
# Запуск: PYTHONIOENCODING=utf-8 python make_detector_card.py  → пишет detector_card.CARD_PATH (перезапись — с бэкапом .bak).
import os, sys, json, shutil
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import detector_card as dc

SL = [[238.632, 24.073], [583.187, 40.624], [911.204, 57.351], [964.766, 59.686], [968.971, 59.867], [2614.511, 107.420]]
card = {"detector": "Atom GS2020 PRO MAX", "vessel": "Маринелли 1 л (РАДЭК)", "version": "2026-10-09",
        "reference": {"spectrum": "Th-232 ОИСН-10 в Маринелли 1 л, 1161 мл по уровню, ρ 0,9061 г/см³ (референс 26.09.2026)",
                      "fwhm_source": "СпектраЛайн «Параметры пиков» 2026-09-26 (оператор: «точнее сделал»), скрин в Референсы\\",
                      "scale_source": "#SHAPE-1 26.09: минимум невязки формы пиков референса, множитель 1,05 на всех точках"},
        "fwhm": {"points_sl": SL, "scale": {str(E): 1.05 for E, _ in SL}, "interp": "кусочно-линейная в log(E)–log(ПШПВ) (mix_unfold_g1s.make_fwhm)",
                 "excluded": [{"E": 1460.822, "fwhm_sl": 70.963, "why": "пик K-40 фона в спектре тория — не свойство референсной пробы; в свёртку не идёт"}]},
        "peak_shape": {"kernel": "gauss", "tail_T": 0.0, "blur": 1.0, "why": "хвост донора Гамма-1С (T=0,75) сдвигал пики до +4,8 кэВ (W-149)"},
        "not_in_card": "шкала энергии, живое время, фон — свойства каждого замера (#CAL-0)"}
card["fwhm"]["sha256"] = dc.points_sha(dc.effective_points(card))
os.makedirs(os.path.dirname(dc.CARD_PATH), exist_ok=True)
if os.path.exists(dc.CARD_PATH): shutil.copy2(dc.CARD_PATH, dc.CARD_PATH + ".bak")
json.dump(card, open(dc.CARD_PATH, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("карточка:", dc.CARD_PATH); print("точки свёртки:", dc.effective_points(card)); print("sha256:", card["fwhm"]["sha256"])
