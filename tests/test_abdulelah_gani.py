from pathlib import Path

import numpy as np

import pandas as pd

import pytest

from ugropy import abdulelah_gani


# =============================================================================
# Util function
# =============================================================================
def pretty_sol(groups_vector: np.ndarray):
    p_idx = np.linspace(0, 219, 220, dtype=int)
    s_idx = np.linspace(220, 349, 130, dtype=int)
    t_idx = np.linspace(350, 423, 74, dtype=int)

    p_occurs = groups_vector[p_idx]
    s_ocurrs = groups_vector[s_idx]
    t_ocurrs = groups_vector[t_idx]

    p_dict = {}
    dfp = abdulelah_gani.primary_model.subgroups_info
    for p, pocc in zip(p_idx, p_occurs):
        if pocc > 0:
            group_name = dfp[dfp["group_number"] == p + 1].index[0]
            p_dict[group_name] = int(pocc)

    s_dict = {}
    dfs = abdulelah_gani.secondary_model.subgroups_info
    for s, s_occ in zip(s_idx, s_ocurrs):
        if s_occ > 0:
            group_name = dfs[dfs["group_number"] == s + 1].index[0]
            s_dict[group_name] = int(s_occ)

    t_dict = {}
    dft = abdulelah_gani.tertiary_model.subgroups_info
    for t, t_occ in zip(t_idx, t_ocurrs):
        if t_occ > 0:
            group_name = dft[dft["group_number"] == t + 1].index[0]
            t_dict[group_name] = int(t_occ)

    return p_dict, s_dict, t_dict


# =============================================================================
# Data path and groups vector
# =============================================================================
_here = here = Path(__file__).parent
dfs_dir = _here / "abdulelah_gani_frags"

groups = np.linspace(1, 424, 424, dtype=int).astype(str)

# =============================================================================
# Critical temperature
# =============================================================================
df_tc = pd.read_csv(
    dfs_dir / "tc.csv", index_col="SMILES", sep="|", comment="?"
)
df_tc.dropna(inplace=True)


@pytest.mark.agani
@pytest.mark.parametrize("smiles", df_tc.index)
def test_abdulelah_gani_tc(smiles):
    sols = abdulelah_gani.get_groups(
        smiles, "smiles", search_multiple_solutions=True
    )

    ag_sol = df_tc.loc[smiles, groups].values.astype(int)

    if not any([np.allclose(sol.ml_vector.ravel(), ag_sol) for sol in sols]):
        # Expected solutions
        ep, es, et = pretty_sol(ag_sol)

        # Obtained solutions
        op = [s.primary.subgroups for s in sols]
        obtained_primary = "\n".join(str(o) for o in op)

        error_message = (
            "Critical Temperature Abdulelah-Gani model failed for\n"
            f"SMILES: {smiles}\n"
            f"Expected primary:\n {str(ep)}\n"
            f"Expected secondary:\n {str(es)}\n"
            f"Expected tertiary:\n {str(et)}\n"
            "------------------------------------\n"
            "obtained primary:\n" + obtained_primary + "\n"
            f"obtained secondary:\n {str(sols[0].secondary.subgroups)}\n"
            f"obtained tertiary:\n {str(sols[0].tertiary.subgroups)}\n"
        )

        assert False, error_message

    for sol in sols:
        if np.allclose(sol.ml_vector.ravel(), ag_sol):
            agani_tc = df_tc.loc[smiles, "GC-SIMPLE"]
            assert np.isclose(
                sol.critical_temperature.to("K").magnitude, agani_tc
            ), (
                f"Critical Temperature Abdulelah-Gani model failed for\n"
                f"SMILES: {smiles}\n"
                f"Expected: {agani_tc}\n"
                f"Obtained: {sol.critical_temperature.magnitude}"
            )


# =============================================================================
# Critical pressure
# =============================================================================
df_pc = pd.read_csv(
    dfs_dir / "pc.csv", index_col="SMILES", sep="|", comment="?"
)
df_pc.dropna(inplace=True)


@pytest.mark.agani
@pytest.mark.parametrize("smiles", df_pc.index)
def test_abdulelah_gani_pc(smiles):
    sols = abdulelah_gani.get_groups(
        smiles, "smiles", search_multiple_solutions=True
    )

    ag_sol = df_pc.loc[smiles, groups].values.astype(int)

    if not any([np.allclose(sol.ml_vector.ravel(), ag_sol) for sol in sols]):
        # Expected solutions
        ep, es, et = pretty_sol(ag_sol)

        # Obtained solutions
        op = [s.primary.subgroups for s in sols]
        obtained_primary = "\n".join(str(o) for o in op)

        error_message = (
            "Critical pressure Abdulelah-Gani model failed for\n"
            f"SMILES: {smiles}\n"
            f"Expected primary:\n {str(ep)}\n"
            f"Expected secondary:\n {str(es)}\n"
            f"Expected tertiary:\n {str(et)}\n"
            "------------------------------------\n"
            "obtained primary:\n" + obtained_primary + "\n"
            f"obtained secondary:\n {str(sols[0].secondary.subgroups)}\n"
            f"obtained tertiary:\n {str(sols[0].tertiary.subgroups)}\n"
        )

        assert False, error_message

    for sol in sols:
        if np.allclose(sol.ml_vector.ravel(), ag_sol):
            agani_pc = df_pc.loc[smiles, "GC-Simple"]
            assert np.isclose(
                sol.critical_pressure.to("bar").magnitude, agani_pc
            ), (
                f"Critical pressure Abdulelah-Gani model failed for\n"
                f"SMILES: {smiles}\n"
                f"Expected: {agani_pc}\n"
                f"Obtained: {sol.critical_pressure.magnitude}"
            )


# =============================================================================
# Critical volume
# =============================================================================
df_vc = pd.read_csv(
    dfs_dir / "vc.csv", index_col="SMILES", sep="|", comment="?"
)
df_vc.dropna(inplace=True)


@pytest.mark.agani
@pytest.mark.parametrize("smiles", df_vc.index)
def test_abdulelah_gani_vc(smiles):
    sols = abdulelah_gani.get_groups(
        smiles, "smiles", search_multiple_solutions=True
    )

    ag_sol = df_vc.loc[smiles, groups].values.astype(int)

    if not any([np.allclose(sol.ml_vector.ravel(), ag_sol) for sol in sols]):
        # Expected solutions
        ep, es, et = pretty_sol(ag_sol)

        # Obtained solutions
        op = [s.primary.subgroups for s in sols]
        obtained_primary = "\n".join(str(o) for o in op)

        error_message = (
            "Critical volume Abdulelah-Gani model failed for\n"
            f"SMILES: {smiles}\n"
            f"Expected primary:\n {str(ep)}\n"
            f"Expected secondary:\n {str(es)}\n"
            f"Expected tertiary:\n {str(et)}\n"
            "------------------------------------\n"
            "obtained primary:\n" + obtained_primary + "\n"
            f"obtained secondary:\n {str(sols[0].secondary.subgroups)}\n"
            f"obtained tertiary:\n {str(sols[0].tertiary.subgroups)}\n"
        )

        assert False, error_message

    for sol in sols:
        if np.allclose(sol.ml_vector.ravel(), ag_sol):
            agani_vc = df_vc.loc[smiles, "GC-Simple Prediction"]
            assert np.isclose(
                sol.critical_volume.to("cm^3/mol").magnitude, agani_vc
            ), (
                f"Critical volume Abdulelah-Gani model failed for\n"
                f"SMILES: {smiles}\n"
                f"Expected: {agani_vc}\n"
                f"Obtained: {sol.critical_volume.magnitude}"
            )
            
            
# =============================================================================
# Acentric factor
# =============================================================================
df_w = pd.read_csv(
    dfs_dir / "w.csv", index_col="SMILES", sep="|", comment="?"
)
df_w.dropna(inplace=True)


@pytest.mark.agani
@pytest.mark.parametrize("smiles", df_w.index)
def test_abdulelah_gani_w(smiles):
    sols = abdulelah_gani.get_groups(
        smiles, "smiles", search_multiple_solutions=True
    )

    ag_sol = df_w.loc[smiles, groups].values.astype(int)

    if not any([np.allclose(sol.ml_vector.ravel(), ag_sol) for sol in sols]):
        # Expected solutions
        ep, es, et = pretty_sol(ag_sol)

        # Obtained solutions
        op = [s.primary.subgroups for s in sols]
        obtained_primary = "\n".join(str(o) for o in op)

        error_message = (
            "Acentric factor Abdulelah-Gani model failed for\n"
            f"SMILES: {smiles}\n"
            f"Expected primary:\n {str(ep)}\n"
            f"Expected secondary:\n {str(es)}\n"
            f"Expected tertiary:\n {str(et)}\n"
            "------------------------------------\n"
            "obtained primary:\n" + obtained_primary + "\n"
            f"obtained secondary:\n {str(sols[0].secondary.subgroups)}\n"
            f"obtained tertiary:\n {str(sols[0].tertiary.subgroups)}\n"
        )

        assert False, error_message

    # TODO: Why so unprecise? Some parameters are wrong?
    # for sol in sols:
    #     if np.allclose(sol.ml_vector.ravel(), ag_sol):
    #         agani_w = df_w.loc[smiles, "GC-Simple Prediction"]
    #         assert np.isclose(
    #             sol.acentric_factor.magnitude, agani_w, atol=1e-1
    #         ), (
    #             f"Acentric factor Abdulelah-Gani model failed for\n"
    #             f"SMILES: {smiles}\n"
    #             f"Expected: {agani_w}\n"
    #             f"Obtained: {sol.acentric_factor.magnitude}"
    #         )