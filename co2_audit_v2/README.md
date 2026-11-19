# Krátkovlnná absorpce CO₂: numerický výzkum a audit, rozšířená verze

Hlavní zpráva je `index.html`. Všechny její grafy jsou vložené přímo do HTML,
takže ke čtení nepotřebuje síť. Schválený český úvod je obsahově zachovaný;
na výslovné přání uživatele byly v jeho nadpisu doplněny uvozovky kolem
výrazu „zadržování tepla“. Záhlaví zůstalo beze změny. Původní kopie úvodu
je v `audit_evidence/intro_before_quotation_update.html`, původní i aktuální
kontrolní otisky jsou v `preservation.json`.

Hned za úvodem je nový popularizační souhrn `popular_summary.html`. Vysvětluje
také odlišení radiačního příspěvku od teplotní zpětné vazby a vztah obnovy
ozonu po Montrealském protokolu k teplotním trendům stratosféry. Nové odvozené
teplotní podíly počítá `discussion_metrics.py` z již uložených výsledků;
nejde o nový ozonový nebo klimatický experiment.

Integritu celé dodávky lze před instalací vědeckého prostředí ověřit obyčejným
Pythonem: `python3 verify_manifest.py`. Tento krok nic nestahuje ani nemění.

## Co bylo provedeno

- RTE-RRTMGP: skutečný SW i LW přenos pro 100 vážených profilů RFMIP,
  všech pět požadovaných přechodů koncentrace a tři kontrolní definice.
- FDH řešení teploty stratosféry při fixované troposféře, povrchu a složení.
  Dvojitý rozdíl s vypnutím pouze SW CO₂ odděluje SW vyvolanou LW odezvu.
- Samostatné běhy RRTMG a zkompilovaného SOCRATES6/260 na skutečných
  CKDMIP vstupech a porovnání s nezávislými LBL referenčními toky.
- Globální SW srovnání tří rodin modelů na stejných profilech a denní geometrii.
- Šest AFGL atmosfér a tři další vlhkostní varianty, 0,5km síť do80km;
  úplný přepočet plynných optických vlastností při změně vlhkosti.
- Nízký kapalný a vysoký ledový oblak, předepsaný aerosol, tři albeda;
  u obou oblaků také skutečná FDH adaptace.
- Koncentrace140–2240ppm, oddělené numerické a modelové citlivosti,
  vertikální/spektrální/sluneční konvergence a tropopauza±1hladina.
- Vlastní statistický přepočet čtyř veřejných Allenových metanových simulací,
  40ročních průměrů, autokorelace a blokového bootstrapu.
- Kontrola vybraných autorských Cohenových stavů; audit Shine tabulky3.
- Jedenáct věcných grafů a CSV včetně úplných svislých toků a ohřevů.

Čísla v HTML jsou z uložených výsledků, nikoli ručně zadané výsledky článků.
Článek nebo referenční tok představuje testovaná data, ne podmínku pro úpravu
našeho výstupu. Nebyla aplikována žádná procentní korekce teploty podle chyby RRTMG.

## Hlavní výsledky a jejich význam

Základní RRTMGP224/256 + TSIS pásmové váhy,280→420ppm, přibližný roční
globální **jasný** průměr:

- SW absorbovaný povrchem při pevné teplotě: −0,1708433W/m².
- SW absorpce celé atmosféry: +0,2034214W/m².
- Troposféra: +0,0245122; tropopauza až1hPa:+0,1699353;
  nad1hPa:+0,0089738W/m².
- Po FDH, pouze SW vyvolaná část: LW ven naTOA:+0,1072544W/m²;
  LW dolů na tropopauze:+0,0671423W/m²; čistý tok naTOA:−0,0743302W/m².
- Celkové působení SW+LW po FDH:+2,6116463W/m².

Podrobné znaménkové konvence a kontrolní experimenty jsou v HTML. Jedná se
o radiační příspěvky. Z těchto čísel není odhadováno ochlazení mezní vrstvy v°C.
Rozdíl RRTMGP/SOCRATES pro okamžitý jasný povrchový účinek je přibližně
−0,171 až−0,209W/m², nikoli95% statistický interval.

## Instalace na macOS/Apple Silicon

Soubor `environment-osx-arm64.lock.txt` obsahuje přesné veřejné conda balíky
pro platformu osx-arm64. Použitý Python je3.11.18. Lze použít micromamba:

```sh
micromamba create -y -p ./work/env -f environment-osx-arm64.lock.txt
```

Po aktivaci prostředí nastavte pro přehlednost místní pracovní adresáře:

```sh
export CO2_V2_WORK="$PWD/work/runtime"
export MPLCONFIGDIR="$CO2_V2_WORK/mpl"
export NUMBA_CACHE_DIR="$CO2_V2_WORK/numba"
export PYTHONPYCACHEPREFIX="$CO2_V2_WORK/pycache"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
python acquire_data.py --check
python verify.py
```

`verify.py` kontroluje již uložené výsledky a analyticky ověřuje přenos.
`acquire_data.py --check` je plně offline. Bez `--check` stahuje jen chybějící
veřejné vstupy a odmítne soubor s odlišným SHA256 nebo velikostí.
Změněný živý zdroj musí být zpracován jako nová verze experimentu.

Úplný přepočet:

```sh
python run_all.py
```

Přepočet hlavní RRTMGP sady bez nové kompilace SOCRATES:

```sh
python run_all.py --core
```

Uložené NetCDF výsledky FDH se znovu používají. Chcete-li skutečně znovu
počítat konkrétní experiment, odstraňte nebo přejmenujte jeho adresář výsledků.
Mějte volné pracovní místo pro Fortranové zdroje a kompilaci. Neobsahuje to
žádnou novou CESM2 integraci. Na jiné platformě vytvořte odpovídající prostředí
ze seznamu verzí v `environment.json`; binární osx-arm64 lock není přenositelný.

### SOCRATES

`setup_socrates.py` použije připnutý zdrojový archiv dodaný v `vendor/`
(při jeho absenci stáhne veřejnou kopii), kontroluje jeho
otisk, nastaví GNU Fortran a cesty NetCDF v `make/Mk_cmd` a spustí
`build_code l_run_cdf`. Pro malé upstream C utility potřebuje Apple command-line
tools. Zdrojový snapshot je48a7a19; přesný archiv a SHA256 jsou v
`model_sources.json`. Jeho BSD copyright zůstává v `audit_evidence` i zdroji.
Přibalený je také malý původní zdrojový archiv pyRTE-RRTMGP0.2.5; instalované
binární vědecké závislosti jsou přesně připnuté v conda locku.

`SOCRATES_DIR` může ukazovat na již sestavený strom. `socrates_benchmark.py`
používá upstream nctools a Cl_run_cdf, ne emulaci radiačního solveru.
Do pracovních adresářů ukládá vstupy a log pro každý úhel a koncentraci.
Původní SW regresní toky byly reprodukovány s maximální odchylkou0,000122W/m².
LW odchylka proti starším KGO až0,184W/m² je zveřejněna; hlavní LW-FDH
výsledky proto pocházejí z RRTMGP.

## Skripty a výstupy

| Soubor | Úloha |
|---|---|
| `radiation.py` | RFMIP převody jednotek, geometrie, výšky, tropopauza, RTE-RRTMGP, FDH |
| `rrtmg_adapter.py` | Přímé volání AER Fortran SW/LW kernelů |
| `benchmark_ckdmip.py` | Nové RRTMG/RRTMGP výpočty proti skutečným LBL referencím |
| `socrates_benchmark.py` | SOCRATES6/260 na stejných CKDMIP sloupcích |
| `global_socrates.py` | Stejná globální denní geometrie jako RRTMGP |
| `global_rrtmg.py` | Výslovně méně přesná FDH kontrola s tolerancí0,02K/den |
| `scenarios.py` | AFGL, vlhkost, oblaky, aerosol, albedo, koncentrace |
| `fdh_sensitivities.py` | Skutečná FDH pro dvě tropopauzy a dva oblaky |
| `audit_allen_data.py` | Plošné a časové váhy, Welch a blokový bootstrap |
| `check_cohen.py` | Kontrola vybraných uložených autorských stavů |
| `verify.py` | Analytické kontroly, energie, jednotky a rezidua |
| `analyze.py` | CSV, 11grafů, modelové chyby výškového ohřevu |
| `discussion_metrics.py` | Odvozené teplotní podíly a srovnání s5,35ln ze zachovaných běhů |
| `popular_summary.html` | Popularizační text s vysvětlením zavádějících zkratek |
| `build_report.py` | Samostatné HTML s úvodem a požadovanou úpravou uvozovek |
| `data_sources.json` | URL, velikost a SHA256 všech dodaných hlavních vstupů |
| `experiments.json` | Přehled konfigurací a jejich odpovídajících skriptů |

### Datové konvence

Všechny tlaky hlavního modelu jsou vPa, teploty vK, toky veW/m² na
vodorovnou plochu. Pozitivní čistý tok je dolů minus nahoru. Pojmenované
jednosměrné toky jsou kladné ve svém směru. Vrstvy a rozhraní jsou řazeny
od vrcholu k povrchu. Kladná absorpce vrstvy je Fnet(horní)−Fnet(dolní).
K/den je tato energie dělená hmotnostním tepelným obsahem vrstvy.

CSV `global_height_profiles.csv` obsahuje identifikaci profilu, váhu, tlak,
výšku, okamžitý i adaptovaný SW/LW příkon a K/den pro obě kontrolní rodiny.
Pro SW vyvolaný profil je třeba od sebe odečíst řádky `full` a `co2off`.
`absolute_flux_profiles.csv` obsahuje absolutní up/down/direct toky pro čtyři
hlavní koncentrace. Kompletní NetCDF obsahuje i ostatní požadované koncentrace.
`model_heating_profiles.csv` obsahuje všechny modely,50profilů,5úhlů a vrstvy
CKDMIP; statistika se hodnotí jen od tlaku2Pa dolů podle protokolu.

NetCDF otevíraný s automatickým dekódováním převádí čas RFMIP na datum.
Pro surové vstupní číslo používejte `decode_cf=False`; geometrie podporuje
i datum a převádí je na den roku. V xarray se zde záměrně nepřidávají
neslučitelné indexové souřadnice `layer` a `level`, protože by při diferenčním
výpočtu hmotnosti mohly posunout vrstvy při automatickém zarovnání.

## Limity a nalezené problémy

1. Nový úplný globální LBL z HITRAN2024 není dokončenou referencí. Hlavní
   metoda je otevřený CKD model s nezávislým LBL ověřením. Koeficienty datv1.9
   mají metadata2018/2019 a nejsou nově vytvořené z databáze2024.
2. V RRTMGP224 byl nalezen podezřelý kladný povrchový příspěvek pásma
   3,73–12,20µm:+0,00738W/m². Je zachován a výslovně označen v HTML
   i diagnostice. Energetická uzávěrka neznamená spektrální správnost.
3. Horní tlak globálních profilů byl omezen asi na1Pa; modelová hranice je
   podle sloupce73–82km. Nejde o validovaný non-LTE model nad80km.
4. Zjemnění původních60vrstev neobnovuje nezachycenou vertikální variabilitu.
   RFMIP váhy a naše denní kvadratura dávají přibližný jasný roční průměr,
   nikoli nezávislý náhodný vzorek pro výpočet95% globálního intervalu.
5. Přísná FDH konvergence RRTMG se nezdařila; log je zachován. Dodatečný
   méně přesný výsledek má vlastní uvedenou toleranci. RRTMGP hlavní sada
   má rezidua pod2×10⁻⁵K/den a přímou kontrolu obou hranic.
6. Oblačnost, aerosol a albeda jsou kontrolované předepsané scénáře. Globální
   skutečné all-sky váhy ani oblačné/dynamické rychlé úpravy nejsou odvozené.
7. Zenodo10357888 neobsahuje CO₂ běhy. Přiložená metadata to dokazují.
   Vlastní přepočet CH₄ dat není vydáván za reprodukci CO₂ −0,38K.
8. Statistický interval Allenových dat a modelový rozptyl našich toků jsou
   různé druhy nejistoty. Nesměšují se.

Původní omezený čárový experiment je zabalen samostatně v
`legacy/co2_audit_v1.zip`, aby byl reprodukovatelný zachovaný příklad úvodu.
Jeho report popisuje starší rozsah; aktuální závěry jsou v této verzi.

## Zdroje

Primární zdroje a DOI jsou uvedeny přímo u tvrzení v HTML. Datové URL a otisky
jsou v `data_sources.json`, připnuté zdrojové archivy v `model_sources.json`,
přesné balíky v conda locku. Zahrnuté modely a data mají své upstream licence;
žádné členství ani přístupové údaje nebyly součástí balíku. Měření ani
publikované klimatické simulace nejsou označeny za vlastní měření přírody.
