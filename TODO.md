# To do: atom library (open items from the Rietimat / Rietipy work)

Written 2026-10-21 by the Rietimat session. Each item was tested against the current code. Rietimat's checks for the
library are in `Autofit_spectra/Rietimat/test_atom_library_fixes.m` (run it after changing anything here).

## Done in this round
- Label rule (decided by the owner): every form of `Hs`, `He`, `Ho` and `Hw` is hydrogen (`Ho3+`, `He1`, `hs` too), `Oc` is neutral
  oxygen, `Nc` is nitrogen. `Hf` (hafnium) and `Hg` (mercury) are the real elements. In `atomic_scattering_factors.m` the H rule
  now excludes `hf` and `hg` (done in `Github_atom/atom` and `functions/atom`). `element_atom.m` already did this (`He`, `Ho`, `Hs`
  to H, `Hw` kept as `water_H`, `Hf` kept, `Oc`/`Nc` to O/N), so no change was needed there.

## Open
0. **A bare force-field aluminium label is neutral Al.** `atomic_scattering_factors('Alo')` (and any `Al*` other than exactly `Al` and `Alt`)
   falls through the table lookup, prints "Did not found the atomtype", cuts the label to `Al` and uses neutral aluminium (13 electrons).
   `Al` and `Alt` give `Al1.5+` (12.5), as does Rietipy for `Alo`. Through `xrd_atom` this does not matter because `element_atom` turns `Alo`
   into `Al` first. Seen in Rietimat's new scattering table (`res.scattering`) when a phase carries `Alo` labels directly.
1. (done 2026-10-21: `element_atom.m` now keeps `Hg` as mercury, in `Github_atom/atom` and `functions/atom`.)
2. **`atom_MATLAB_toolbox/general_functions/atomic_scattering_factors.m` inside the atomipy repository is stale**
   (`Coding/Windsurf/atominpython/atom_MATLAB_toolbox/general_functions/`). It is an earlier patch: protects only He/Hf, replaces
   `O2-`, is case sensitive. Copy the canonical file over it. The other three patched files there match the canonical ones.
3. **Old copies keep the old rules.** Of 1,308 copies of `atomic_scattering_factors.m` under `MATLAB/` and `Coding/`, 60 are
   identical to the pre-fix file (`functions/atom_3.0`, `Coding/Windsurf/XRD_IN_ATOM`, `Coding/Antigravity/atominpython/atom_MATLAB_functions`,
   many Kebnekaise job folders). `functions/atom_2.12` and `Matlab_drive/atom` differ in all four patched files. Decide which to retire.
4. **`xrd_atom_legacy.m`** was not touched: it has the same file-name defect as `xrd_atom('file.gro')` had (uses `Box_dim` before it
   is defined) and always writes `xrd.dat`.
5. **Unknown labels** still print and fall back by first letter. A `warning` with an identifier would be better (atomipy warns).
6. **No MATLAB counterpart of atomipy's `condense_supercell`** (sites, occupancies and B from the replicas of an MD box).
   `rietveld_structure.m` in Rietimat de-replicates by folding; a library function would be useful.
7. Optional: anomalous dispersion (f', f'') is not in the MATLAB library; Rietipy has it (gemmi Cromer-Liberman). A tabulated
   version would let `xrd_atom` do Fe/Co/Cu-radiation work properly.
