from cases import TCase

import pytest

from ugropy import psrk


@pytest.mark.dependency(name="psrk")
@pytest.mark.psrk
class TestPSRK(TCase):
    # Store al the groups detected in all the cases here:
    tested_groups = set()

    def asserts(self, case, solver):
        if case.psrk_result is None:
            pytest.skip(
                f"No PSRK result defined for {case.identifier}, "
                f"{case.cases_module}"
            )

        result = psrk.get_groups(
            identifier=case.identifier,
            identifier_type=case.identifier_type,
            solver=solver,
            search_multiple_solutions=True,
        )

        if len(result) > 1:
            for r in result:
                comp = r.subgroups in case.psrk_result

                if not comp:
                    message = (
                        "\n"
                        f"Case: {case.identifier}\n"
                        f"Test module: {case.cases_module}\n"
                        f"Expected: {case.psrk_result}\n"
                        f"Obtained: {[r.subgroups for r in result]}"
                    )

                    assert False, message

                self.tested_groups.update(r.subgroups.keys())
        else:
            comp = result[0].subgroups == case.psrk_result

            self.tested_groups.update(result[0].subgroups.keys())

            if not comp:
                message = (
                    "\n"
                    f"Case: {case.identifier}\n"
                    f"Test module: {case.cases_module}\n"
                    f"Expected: {case.psrk_result}\n"
                    f"Obtained: {result[0].subgroups}"
                )

                assert False, message


@pytest.mark.dependency(depends=["psrk"])
@pytest.mark.psrk
def test_psrk_groups_coverage():
    # Check if all the groups were detected on at least one case
    for group in TestPSRK.tested_groups:
        assert group in psrk.subgroups.index, f"Group {group} not tested"
