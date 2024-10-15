from cases import (
    acids_cases,
    alcohols_cases,
    aldehydes_cases,
    aromatics_cases,
    Case,
    TCase,
    complex_cases,
    epoxides_cases,
    esters_cases,
    ethers_cases,
    halogens_cases,
    hydrocarbons_cases,
    insaturated_hydrocarbons_cases,
    ketones_cases,
    nitrogen_cases,
    particulars_cases,
    silicon_cases,
    sulfur_cases,
)

import pytest

from ugropy import unifac


class TestUNIFAC:
    @pytest.mark.parametrize("case", hydrocarbons_cases)
    def test_unifac_hydrocarbons(self, case: Case):
        if case.unifac_result is None:
            pytest.skip(
                f"No UNIFAC result for {case.identifier}, {case.cases_module}"
            )
        else:
            result = unifac.get_groups(
                case.identifier,
                case.identifier_type,
                search_multiple_solutions=True,
            )

            if len(result) > 1:
                assert (
                    False
                ), f"Fijate este tiene multiples: {case.identifier}, {case.cases_module}"

            same_groups = result[0].subgroups == case.unifac_result
